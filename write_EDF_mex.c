/**
 * WRITE_EDF_MEX  High-performance EDF / EDF+ writer for MATLAB
 *
 * Mirror of read_EDF_mex.c. Streams records one at a time and supports
 * gzip output (.edf.gz) on the fly via vendored zlib.
 *
 * MATLAB usage:
 *   write_EDF_mex(filename, header, signal_header, signal_cell, annotations,
 *                 autoscale_mode, verbose, debug)
 *
 *     autoscale_mode : 0 = preserve (use existing physical_min/max, clip)
 *                      1 = recompute (set physical_min/max from data)
 *
 *   The MATLAB side validates inputs before calling. Channels and signal
 *   headers must be consistent (every non-annotation signal must have
 *   length == samples_in_record * num_data_records).
 *
 * Compilation:
 *   mex -O -largeArrayDims -Izlib write_EDF_mex.c zlib/(asterisk).c
 */

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>

#if defined(MX_COMPAT_32)
#error "This MEX requires -largeArrayDims"
#endif

#define EDF_HEADER_SIZE 256

static int g_debug = 0;
#define DBG(...) do { if (g_debug) mexPrintf(__VA_ARGS__); } while (0)

/* ============================================================
 *                      WRITER IO ABSTRACTION
 * ============================================================ */

typedef struct {
    int     is_gz;
    FILE   *fp;
    gzFile  gz;
} EDFWriter;

static int has_gz_suffix(const char *fname)
{
    size_t n = strlen(fname);
    if (n < 3) return 0;
    return (tolower((unsigned char)fname[n-3]) == '.' &&
            tolower((unsigned char)fname[n-2]) == 'g' &&
            tolower((unsigned char)fname[n-1]) == 'z');
}

static int edfw_open(EDFWriter *w, const char *fname)
{
    w->fp = NULL;
    w->gz = NULL;
    w->is_gz = has_gz_suffix(fname);

    if (w->is_gz) {
        w->gz = gzopen(fname, "wb");
        return (w->gz != NULL);
    } else {
        w->fp = fopen(fname, "wb");
        return (w->fp != NULL);
    }
}

static int edfw_write(EDFWriter *w, const void *buf, size_t n)
{
    if (w->is_gz) {
        int put = gzwrite(w->gz, buf, (unsigned)n);
        return (put == (int)n);
    } else {
        return (fwrite(buf, 1, n, w->fp) == n);
    }
}

static void edfw_close(EDFWriter *w)
{
    if (w->is_gz) {
        if (w->gz) { gzclose(w->gz); w->gz = NULL; }
    } else {
        if (w->fp) { fclose(w->fp); w->fp = NULL; }
    }
}

/* ============================================================
 *                    HEADER FIELD HELPERS
 * ============================================================
 *
 * EDF text fields are space-padded ASCII, left-justified, fixed width.
 * Truncate over-long content to fit; pad with spaces.
 */

static void pad_field(char *dst, mwSize width, const char *src)
{
    mwSize n = src ? (mwSize)strlen(src) : 0;
    if (n > width) n = width;
    if (n) memcpy(dst, src, n);
    if (n < width) memset(dst + n, ' ', width - n);
}

/* Format an integer into a fixed-width left-justified ASCII field. */
static void pad_int(char *dst, mwSize width, long long v)
{
    char buf[32];
    int n = snprintf(buf, sizeof(buf), "%lld", v);
    if (n < 0) { memset(dst, ' ', width); return; }
    if ((mwSize)n > width) n = (int)width;  /* truncate digits if absurd */
    memcpy(dst, buf, n);
    if ((mwSize)n < width) memset(dst + n, ' ', width - n);
}

/* Format a double into a fixed-width left-justified ASCII field.
 * Strategy: try increasing precision until the value fits without
 * losing more than necessary, then strip trailing zeros after a '.'. */
static void pad_double(char *dst, mwSize width, double v)
{
    char buf[64];
    int n;
    /* Special handling: integers print without decimals if they fit. */
    if (v == floor(v) && fabs(v) < 1e15) {
        long long lv = (long long)v;
        n = snprintf(buf, sizeof(buf), "%lld", lv);
    } else {
        n = snprintf(buf, sizeof(buf), "%.6f", v);
        /* Strip trailing zeros */
        if (n > 0 && strchr(buf, '.')) {
            while (n > 1 && buf[n-1] == '0') n--;
            if (n > 1 && buf[n-1] == '.') n--;
            buf[n] = '\0';
        }
    }
    if (n < 0) { memset(dst, ' ', width); return; }
    if ((mwSize)n > width) n = (int)width;
    memcpy(dst, buf, n);
    if ((mwSize)n < width) memset(dst + n, ' ', width - n);
}

/* ============================================================
 *                STRUCT FIELD READ HELPERS (mxArray)
 * ============================================================ */

static const char *get_string_field(const mxArray *s, mwIndex idx, const char *fname,
                                    char *buf, size_t buflen)
{
    mxArray *f = mxGetField(s, idx, fname);
    if (!f || !mxIsChar(f)) {
        if (buflen) buf[0] = '\0';
        return buf;
    }
    mxGetString(f, buf, (mwSize)buflen);
    return buf;
}

static double get_scalar_field(const mxArray *s, mwIndex idx, const char *fname, double dflt)
{
    mxArray *f = mxGetField(s, idx, fname);
    if (!f || mxGetNumberOfElements(f) < 1) return dflt;
    return mxGetScalar(f);
}

/* ============================================================
 *                      ANNOTATION SERIALIZER
 * ============================================================
 *
 * Each record's annotation channel must start with the mandatory EDF+
 * timestamp TAL: "+<rec_offset_seconds>\x14\x14\x00". User-supplied
 * annotations are appended to record 0 by default; if record 0's
 * annotation buffer fills up, spill into subsequent records.
 *
 * This isn't the only valid EDF+ layout, but it's the simplest that
 * round-trips through read_EDF's annotation parser.
 */

typedef struct {
    char  *bytes;     /* annot_bytes_per_record long; zero-init each record */
    mwSize used;      /* bytes filled */
    mwSize cap;       /* annot_bytes_per_record */
} AnnotRecord;

/* Format a double as a +/- decimal, e.g. "+30" or "+1.234".
 * Mirrors what read_EDF expects ('+' or '-' prefix, then digits/dot). */
static int format_onset(char *buf, size_t buflen, double v)
{
    int n;
    if (v == floor(v) && fabs(v) < 1e15) {
        n = snprintf(buf, buflen, "%+lld", (long long)v);
    } else {
        n = snprintf(buf, buflen, "%+.6f", v);
        if (n > 0 && strchr(buf, '.')) {
            int k = n;
            while (k > 1 && buf[k-1] == '0') k--;
            if (k > 1 && buf[k-1] == '.') k--;
            buf[k] = '\0';
            n = k;
        }
    }
    return n;
}

/* Append a TAL "+onset\x14text\x14\x00" or "+onset\x14\x14\x00".
 * Returns 1 on success, 0 if the buffer is full. */
static int append_tal(AnnotRecord *r, double onset, const char *text)
{
    char onset_buf[64];
    int olen = format_onset(onset_buf, sizeof(onset_buf), onset);
    size_t tlen = text ? strlen(text) : 0;
    /* Total bytes: olen + 0x14 + tlen + 0x14 + 0x00 */
    size_t need = olen + 1 + tlen + 1 + 1;
    if (r->used + need > r->cap) return 0;
    memcpy(r->bytes + r->used, onset_buf, olen); r->used += olen;
    r->bytes[r->used++] = 0x14;
    if (tlen) { memcpy(r->bytes + r->used, text, tlen); r->used += tlen; }
    r->bytes[r->used++] = 0x14;
    r->bytes[r->used++] = 0x00;
    return 1;
}

/* ============================================================
 *                         MEX ENTRY
 * ============================================================ */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    (void)nlhs; (void)plhs;

    if (nrhs < 4)
        mexErrMsgIdAndTxt("write_EDF_mex:Input",
            "Usage: write_EDF_mex(filename, header, signal_header, signal_cell, "
            "annotations, autoscale_mode, verbose, debug)");

    char fname[4096];
    mxGetString(prhs[0], fname, sizeof(fname));

    const mxArray *m_hdr   = prhs[1];
    const mxArray *m_shdrs = prhs[2];
    const mxArray *m_data  = prhs[3];
    const mxArray *m_annot = (nrhs > 4) ? prhs[4] : NULL;
    int autoscale_mode = (nrhs > 5) ? (int)mxGetScalar(prhs[5]) : 0;
    int verbose        = (nrhs > 6) ? (int)mxGetScalar(prhs[6]) : 0;
    g_debug            = (nrhs > 7) ? (int)mxGetScalar(prhs[7]) : 0;

    if (!mxIsStruct(m_hdr))
        mexErrMsgIdAndTxt("write_EDF_mex:Input", "header must be a struct");
    if (!mxIsStruct(m_shdrs))
        mexErrMsgIdAndTxt("write_EDF_mex:Input", "signal_header must be a struct array");
    if (!mxIsCell(m_data))
        mexErrMsgIdAndTxt("write_EDF_mex:Input", "signal_cell must be a cell array");

    int num_signals = (int)mxGetNumberOfElements(m_shdrs);
    int n_data_cells = (int)mxGetNumberOfElements(m_data);
    if (num_signals != n_data_cells)
        mexErrMsgIdAndTxt("write_EDF_mex:Mismatch",
            "signal_header has %d entries but signal_cell has %d", num_signals, n_data_cells);
    if (num_signals < 1)
        mexErrMsgIdAndTxt("write_EDF_mex:Input", "no signals to write");

    /* ------------------ Pull header fields ------------------ */
    char edf_ver[16]              = "0";
    char patient_id[128]          = "X X X X";
    char local_rec_id[128]        = "Startdate X X X X";
    char recording_startdate[16]  = "01.01.01";
    char recording_starttime[16]  = "00.00.00";

    get_string_field(m_hdr, 0, "edf_ver",              edf_ver, sizeof(edf_ver));
    get_string_field(m_hdr, 0, "patient_id",            patient_id, sizeof(patient_id));
    get_string_field(m_hdr, 0, "local_rec_id",          local_rec_id, sizeof(local_rec_id));
    get_string_field(m_hdr, 0, "recording_startdate",   recording_startdate, sizeof(recording_startdate));
    get_string_field(m_hdr, 0, "recording_starttime",   recording_starttime, sizeof(recording_starttime));

    double data_record_duration = get_scalar_field(m_hdr, 0, "data_record_duration", 1.0);
    if (!(data_record_duration > 0))
        mexErrMsgIdAndTxt("write_EDF_mex:Input", "data_record_duration must be > 0");

    /* num_header_bytes is ignored on input; we recompute it below. */

    /* ------------------ Per-signal header info ------------------ */
    typedef struct {
        char    label[17];
        char    transducer[81];
        char    phys_dim[9];
        double  phys_min, phys_max;
        double  dig_min,  dig_max;
        char    prefilt[81];
        int     samples_in_record;
        int     is_annotation;
        const double *data;     /* pointer into MATLAB array, NULL for annotation channel */
        mwSize  data_len;
    } SInfo;

    SInfo *si = mxCalloc(num_signals, sizeof(SInfo));

    for (int i = 0; i < num_signals; i++) {
        get_string_field(m_shdrs, i, "signal_labels",      si[i].label,      sizeof(si[i].label));
        get_string_field(m_shdrs, i, "transducer_type",    si[i].transducer, sizeof(si[i].transducer));
        get_string_field(m_shdrs, i, "physical_dimension", si[i].phys_dim,   sizeof(si[i].phys_dim));
        get_string_field(m_shdrs, i, "prefiltering",       si[i].prefilt,    sizeof(si[i].prefilt));

        si[i].phys_min          = get_scalar_field(m_shdrs, i, "physical_min",      0);
        si[i].phys_max          = get_scalar_field(m_shdrs, i, "physical_max",      0);
        si[i].dig_min           = get_scalar_field(m_shdrs, i, "digital_min",      -32768);
        si[i].dig_max           = get_scalar_field(m_shdrs, i, "digital_max",       32767);
        si[i].samples_in_record = (int)get_scalar_field(m_shdrs, i, "samples_in_record", 0);

        si[i].is_annotation = (strcmp(si[i].label, "EDF Annotations") == 0);

        mxArray *cell_i = mxGetCell(m_data, i);
        if (cell_i && mxIsDouble(cell_i)) {
            si[i].data     = mxGetPr(cell_i);
            si[i].data_len = mxGetNumberOfElements(cell_i);
        } else {
            si[i].data     = NULL;
            si[i].data_len = 0;
        }
    }

    /* ------------------ Determine num_data_records ------------------ */
    /* All non-annotation channels must contribute the same number of
     * complete records. If they disagree, the writer truncates to the
     * minimum (and warns). */
    mwSize num_records = 0;
    int    have_count  = 0;

    for (int i = 0; i < num_signals; i++) {
        if (si[i].is_annotation) continue;
        if (si[i].samples_in_record <= 0)
            mexErrMsgIdAndTxt("write_EDF_mex:Input",
                "signal %d ('%s') has samples_in_record <= 0", i, si[i].label);
        mwSize r = si[i].data_len / (mwSize)si[i].samples_in_record;
        if (!have_count) { num_records = r; have_count = 1; }
        else if (r < num_records) num_records = r;
    }
    if (!have_count || num_records == 0)
        mexErrMsgIdAndTxt("write_EDF_mex:NoData",
            "no complete data records to write (need at least one non-annotation channel "
            "with samples_in_record > 0 and data_len >= samples_in_record)");

    if (verbose) mexPrintf("write_EDF_mex: writing %llu records\n",
                           (unsigned long long)num_records);

    /* For the annotation channel, samples_in_record may be 0 (no annotation
     * channel was prepared). Force a sensible default if user wants one: a
     * 60-byte annotation buffer per record (= samples_in_record of 30). */
    int annot_idx = -1;
    for (int i = 0; i < num_signals; i++) {
        if (si[i].is_annotation) {
            annot_idx = i;
            if (si[i].samples_in_record <= 0) si[i].samples_in_record = 30;
            break;
        }
    }

    /* ------------------ Optional: recompute physical range ------------------ */
    if (autoscale_mode == 1) {
        for (int i = 0; i < num_signals; i++) {
            if (si[i].is_annotation || !si[i].data) continue;
            mwSize n = (mwSize)si[i].samples_in_record * num_records;
            double mn = si[i].data[0], mx = si[i].data[0];
            for (mwSize k = 1; k < n; k++) {
                double v = si[i].data[k];
                if (v < mn) mn = v;
                if (v > mx) mx = v;
            }
            if (mn == mx) { mn -= 1; mx += 1; }  /* avoid zero-range */
            si[i].phys_min = mn;
            si[i].phys_max = mx;
        }
    }

    /* ------------------ Compute headers ------------------ */
    int num_header_bytes = EDF_HEADER_SIZE * (1 + num_signals);

    EDFWriter w;
    if (!edfw_open(&w, fname))
        mexErrMsgIdAndTxt("write_EDF_mex:File", "Cannot open '%s' for writing", fname);

    /* ------ Main 256-byte header ------ */
    char mh[EDF_HEADER_SIZE];
    memset(mh, ' ', EDF_HEADER_SIZE);
    pad_field (mh + 0,   8,  edf_ver);
    pad_field (mh + 8,   80, patient_id);
    pad_field (mh + 88,  80, local_rec_id);
    pad_field (mh + 168, 8,  recording_startdate);
    pad_field (mh + 176, 8,  recording_starttime);
    pad_int   (mh + 184, 8,  num_header_bytes);
    /* mh + 192 .. 235 reserved (44 bytes; spaces) */
    pad_int   (mh + 236, 8,  (long long)num_records);
    pad_double(mh + 244, 8,  data_record_duration);
    pad_int   (mh + 252, 4,  num_signals);

    if (!edfw_write(&w, mh, EDF_HEADER_SIZE)) {
        edfw_close(&w);
        mexErrMsgIdAndTxt("write_EDF_mex:Write", "Failed to write main header");
    }

    /* ------ Per-signal header block ------ */
    /* EDF stores signal headers transposed: all labels first, then all
     * transducers, etc. Build into a single buffer of size 256*N. */
    char *sh = mxCalloc(EDF_HEADER_SIZE * num_signals, 1);
    memset(sh, ' ', EDF_HEADER_SIZE * num_signals);

    char *p = sh;
    for (int i = 0; i < num_signals; i++) { pad_field(p, 16, si[i].label);          p += 16; }
    for (int i = 0; i < num_signals; i++) { pad_field(p, 80, si[i].transducer);     p += 80; }
    for (int i = 0; i < num_signals; i++) { pad_field(p, 8,  si[i].phys_dim);       p += 8;  }
    for (int i = 0; i < num_signals; i++) { pad_double(p, 8, si[i].phys_min);       p += 8;  }
    for (int i = 0; i < num_signals; i++) { pad_double(p, 8, si[i].phys_max);       p += 8;  }
    for (int i = 0; i < num_signals; i++) { pad_int   (p, 8, (long long)si[i].dig_min); p += 8; }
    for (int i = 0; i < num_signals; i++) { pad_int   (p, 8, (long long)si[i].dig_max); p += 8; }
    for (int i = 0; i < num_signals; i++) { pad_field (p, 80, si[i].prefilt);       p += 80; }
    for (int i = 0; i < num_signals; i++) { pad_int   (p, 8, (long long)si[i].samples_in_record); p += 8; }
    for (int i = 0; i < num_signals; i++) { /* 32-byte reserved: spaces */          p += 32; }

    if (!edfw_write(&w, sh, EDF_HEADER_SIZE * num_signals)) {
        mxFree(sh);
        edfw_close(&w);
        mexErrMsgIdAndTxt("write_EDF_mex:Write", "Failed to write signal headers");
    }
    mxFree(sh);

    /* ------------------ Per-signal byte offsets within a record ------------------ */
    mwSize total_samp_per_rec = 0;
    for (int i = 0; i < num_signals; i++) total_samp_per_rec += si[i].samples_in_record;
    mwSize bytes_per_rec = total_samp_per_rec * 2;

    mwSize *sig_offset_samples = mxMalloc(num_signals * sizeof(mwSize));
    {
        mwSize acc = 0;
        for (int i = 0; i < num_signals; i++) {
            sig_offset_samples[i] = acc;
            acc += si[i].samples_in_record;
        }
    }

    /* ------------------ Pre-build annotation per-record buffers ------------------ */
    AnnotRecord *arecs = NULL;
    mwSize annot_bytes_per_rec = 0;
    if (annot_idx >= 0) {
        annot_bytes_per_rec = (mwSize)si[annot_idx].samples_in_record * 2;
        arecs = mxCalloc(num_records, sizeof(AnnotRecord));
        for (mwSize r = 0; r < num_records; r++) {
            arecs[r].bytes = mxCalloc(annot_bytes_per_rec, 1);
            arecs[r].cap   = annot_bytes_per_rec;
            arecs[r].used  = 0;
            /* Mandatory record-start timestamp */
            double t = (double)r * data_record_duration;
            char ts_buf[64];
            int olen = format_onset(ts_buf, sizeof(ts_buf), t);
            /* "+t\x14\x14\x00" */
            if ((mwSize)(olen + 3) <= annot_bytes_per_rec) {
                memcpy(arecs[r].bytes, ts_buf, olen);
                arecs[r].bytes[olen]   = 0x14;
                arecs[r].bytes[olen+1] = 0x14;
                arecs[r].bytes[olen+2] = 0x00;
                arecs[r].used = olen + 3;
            }
        }

        /* Distribute user annotations. Place each in the record that matches
         * its onset; spill forward if that record is full. */
        if (m_annot && mxIsStruct(m_annot)) {
            mwSize n_annot = mxGetNumberOfElements(m_annot);
            int dropped = 0;
            for (mwSize a = 0; a < n_annot; a++) {
                double onset = get_scalar_field(m_annot, a, "onset", 0);
                mwSize target = (mwSize)floor(onset / data_record_duration);
                if (target >= num_records) target = num_records - 1;

                /* The 'text' field is a cell array of strings. */
                mxArray *t = mxGetField(m_annot, a, "text");
                int placed = 0;
                if (t && mxIsCell(t)) {
                    mwSize nt = mxGetNumberOfElements(t);
                    for (mwSize k = 0; k < nt; k++) {
                        mxArray *s = mxGetCell(t, k);
                        if (!s || !mxIsChar(s)) continue;
                        char tbuf[256];
                        mxGetString(s, tbuf, sizeof(tbuf));
                        /* Try target, then spill forward. */
                        mwSize r;
                        for (r = target; r < num_records; r++) {
                            if (append_tal(&arecs[r], onset, tbuf)) { placed = 1; break; }
                        }
                        if (!placed) { dropped++; break; }
                    }
                } else if (t && mxIsChar(t)) {
                    char tbuf[256];
                    mxGetString(t, tbuf, sizeof(tbuf));
                    mwSize r;
                    for (r = target; r < num_records; r++) {
                        if (append_tal(&arecs[r], onset, tbuf)) { placed = 1; break; }
                    }
                    if (!placed) dropped++;
                }
            }
            if (dropped && verbose)
                mexPrintf("write_EDF_mex: %d annotation(s) dropped (no record had "
                          "room; consider increasing the annotation channel's "
                          "samples_in_record).\n", dropped);
        }
    }

    /* ------------------ Stream records ------------------ */
    unsigned char *rec_buf = mxCalloc(bytes_per_rec, 1);

    for (mwSize r = 0; r < num_records; r++) {
        memset(rec_buf, 0, bytes_per_rec);

        for (int i = 0; i < num_signals; i++) {
            mwSize spr = (mwSize)si[i].samples_in_record;
            mwSize byte_off = sig_offset_samples[i] * 2;
            unsigned char *dst = rec_buf + byte_off;

            if (si[i].is_annotation) {
                if (annot_idx == i && arecs) {
                    mwSize n = arecs[r].used;
                    if (n > annot_bytes_per_rec) n = annot_bytes_per_rec;
                    memcpy(dst, arecs[r].bytes, n);
                    /* tail already zeroed via memset above */
                }
                continue;
            }

            if (!si[i].data) continue;  /* nothing to write */

            double scale = (si[i].phys_max - si[i].phys_min) /
                           (si[i].dig_max  - si[i].dig_min);
            if (!(scale > 0) && !(scale < 0)) {
                /* zero-range; fall back to ratio of 1 */
                scale = 1.0;
            }
            double inv_scale = 1.0 / scale;
            double dig_min = si[i].dig_min;
            double dig_max = si[i].dig_max;
            double phys_min = si[i].phys_min;

            const double *src = si[i].data + r * spr;
            for (mwSize k = 0; k < spr; k++) {
                double phys = src[k];
                double dig  = (phys - phys_min) * inv_scale + dig_min;
                /* round to nearest, ties away from zero */
                long long iv = (long long)(dig >= 0 ? dig + 0.5 : dig - 0.5);
                if (iv < (long long)dig_min) iv = (long long)dig_min;
                if (iv > (long long)dig_max) iv = (long long)dig_max;
                if (iv < INT16_MIN) iv = INT16_MIN;
                if (iv > INT16_MAX) iv = INT16_MAX;
                int16_t v = (int16_t)iv;
                dst[2*k]   = (unsigned char)(v & 0xFF);
                dst[2*k+1] = (unsigned char)((v >> 8) & 0xFF);
            }
        }

        if (!edfw_write(&w, rec_buf, bytes_per_rec)) {
            mxFree(rec_buf);
            if (arecs) {
                for (mwSize r2 = 0; r2 < num_records; r2++) mxFree(arecs[r2].bytes);
                mxFree(arecs);
            }
            mxFree(sig_offset_samples);
            mxFree(si);
            edfw_close(&w);
            mexErrMsgIdAndTxt("write_EDF_mex:Write",
                              "Failed to write record %llu", (unsigned long long)r);
        }
    }

    mxFree(rec_buf);
    if (arecs) {
        for (mwSize r = 0; r < num_records; r++) mxFree(arecs[r].bytes);
        mxFree(arecs);
    }
    mxFree(sig_offset_samples);
    mxFree(si);

    edfw_close(&w);

    if (verbose) mexPrintf("write_EDF_mex: done (%d signals, %llu records, %s)\n",
                           num_signals, (unsigned long long)num_records,
                           w.is_gz ? ".edf.gz" : ".edf");
}
