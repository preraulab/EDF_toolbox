/*
* =============================================================================
* READ_EDF_MEX.C - Fast MEX File for Reading EDF/EDF+ Files
* =============================================================================
*
* USAGE:
* [header, signalHeader, signalCell, annotations] = read_EDF_mex(edfFN, signalLabels, epochs, verbose, repairHeader)
*
* FEATURES:
* - Full little/big-endian support
* - Auto-repair invalid num_data_records (-1)
* - Robust scaling (handles digital_min == digital_max)
* - EDF+ annotation support (TAL parsing)
*
* COMPILE:
*   mex read_EDF_mex.c   % Linux/macOS/Windows
*
* =============================================================================
*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define EDF_HEADER_SIZE 256
#define MAX_SIGNALS 256

typedef struct {
    char edf_ver[9];
    char patient_id[81];
    char local_rec_id[81];
    char recording_startdate[9];
    char recording_starttime[9];
    int num_header_bytes;
    int num_data_records;
    double data_record_duration;
    int num_signals;
} EDF_Header;

typedef struct {
    char signal_labels[17];
    char transducer_type[81];
    char physical_dimension[9];
    double physical_min;
    double physical_max;
    double digital_min;
    double digital_max;
    char prefiltering[81];
    int samples_in_record;
    double sampling_frequency;
} Signal_Header;

typedef struct {
    double onset;
    mxArray *texts;  /* cell array of strings */
} Annotation;

/* --------------------------------------------------------------- */
/* Utility functions                                               */
/* --------------------------------------------------------------- */
int is_little_endian(void) {
    uint16_t test = 0x0001;
    return *(uint8_t*)&test == 0x01;
}

void trim_string(char *str) {
    int i = 0, j;
    while (str[i] == ' ') i++;
    j = 0;
    while (str[i] != '\0') str[j++] = str[i++];
    while (j > 0 && str[j-1] == ' ') j--;
    str[j] = '\0';
}

int read_edf_header(FILE *fid, EDF_Header *header) {
    unsigned char buffer[EDF_HEADER_SIZE];
    char temp[256];
    if (fread(buffer, 1, EDF_HEADER_SIZE, fid) != EDF_HEADER_SIZE) return 0;

    memcpy(header->edf_ver, buffer + 0, 8);   header->edf_ver[8] = '\0'; trim_string(header->edf_ver);
    memcpy(header->patient_id, buffer + 8, 80); header->patient_id[80] = '\0'; trim_string(header->patient_id);
    memcpy(header->local_rec_id, buffer + 88, 80); header->local_rec_id[80] = '\0'; trim_string(header->local_rec_id);
    memcpy(header->recording_startdate, buffer + 168, 8); header->recording_startdate[8] = '\0'; trim_string(header->recording_startdate);
    memcpy(header->recording_starttime, buffer + 176, 8); header->recording_starttime[8] = '\0'; trim_string(header->recording_starttime);

    memcpy(temp, buffer + 184, 8); temp[8] = '\0'; trim_string(temp); header->num_header_bytes = atoi(temp);
    memcpy(temp, buffer + 236, 8); temp[8] = '\0'; trim_string(temp); header->num_data_records = atoi(temp);
    memcpy(temp, buffer + 244, 8); temp[8] = '\0'; trim_string(temp); header->data_record_duration = atof(temp);
    memcpy(temp, buffer + 252, 4); temp[4] = '\0'; trim_string(temp); header->num_signals = atoi(temp);
    return 1;
}

int read_signal_headers(FILE *fid, Signal_Header *sig_headers, int num_signals) {
    char temp[256];
    int s;

    /* Labels */
    for (s = 0; s < num_signals; s++) {
        unsigned char buf[16];
        if (fread(buf, 1, 16, fid) != 16) return 0;
        memcpy(sig_headers[s].signal_labels, buf, 16);
        sig_headers[s].signal_labels[16] = '\0';
        trim_string(sig_headers[s].signal_labels);
    }
    /* Transducer */
    for (s = 0; s < num_signals; s++) {
        unsigned char buf[80];
        if (fread(buf, 1, 80, fid) != 80) return 0;
        memcpy(sig_headers[s].transducer_type, buf, 80);
        sig_headers[s].transducer_type[80] = '\0';
        trim_string(sig_headers[s].transducer_type);
    }
    /* Physical dimension */
    for (s = 0; s < num_signals; s++) {
        unsigned char buf[8];
        if (fread(buf, 1, 8, fid) != 8) return 0;
        memcpy(sig_headers[s].physical_dimension, buf, 8);
        sig_headers[s].physical_dimension[8] = '\0';
        trim_string(sig_headers[s].physical_dimension);
    }
    /* Physical min/max */
    for (s = 0; s < num_signals; s++) {
        if (fread(temp, 1, 8, fid) != 8) return 0;
        temp[8] = '\0'; trim_string(temp); sig_headers[s].physical_min = atof(temp);
    }
    for (s = 0; s < num_signals; s++) {
        if (fread(temp, 1, 8, fid) != 8) return 0;
        temp[8] = '\0'; trim_string(temp); sig_headers[s].physical_max = atof(temp);
    }
    /* Digital min/max */
    for (s = 0; s < num_signals; s++) {
        if (fread(temp, 1, 8, fid) != 8) return 0;
        temp[8] = '\0'; trim_string(temp); sig_headers[s].digital_min = atof(temp);
    }
    for (s = 0; s < num_signals; s++) {
        if (fread(temp, 1, 8, fid) != 8) return 0;
        temp[8] = '\0'; trim_string(temp); sig_headers[s].digital_max = atof(temp);
    }
    /* Prefiltering */
    for (s = 0; s < num_signals; s++) {
        unsigned char buf[80];
        if (fread(buf, 1, 80, fid) != 80) return 0;
        memcpy(sig_headers[s].prefiltering, buf, 80);
        sig_headers[s].prefiltering[80] = '\0';
        trim_string(sig_headers[s].prefiltering);
    }
    /* Samples per record */
    for (s = 0; s < num_signals; s++) {
        if (fread(temp, 1, 8, fid) != 8) return 0;
        temp[8] = '\0'; trim_string(temp); sig_headers[s].samples_in_record = atoi(temp);
    }
    /* Reserved (32 bytes per signal) */
    for (s = 0; s < num_signals; s++) {
        unsigned char buf[32];
        if (fread(buf, 1, 32, fid) != 32) return 0;
    }
    return 1;
}

int signal_index_from_label(Signal_Header *sig_headers, int num_signals, const char *label) {
    for (int i = 0; i < num_signals; i++) {
        if (strcmp(sig_headers[i].signal_labels, label) == 0) return i;
    }
    return -1;
}

void save_repaired_header(const char *edfFN, int num_records, int verbose) {
    char fixed_filename[4096];
    strcpy(fixed_filename, edfFN);
    char *dot = strrchr(fixed_filename, '.');
    if (dot) {
        char ext[256];
        strcpy(ext, dot);
        *dot = '\0';
        strcat(fixed_filename, "_fixed");
        strcat(fixed_filename, ext);
    } else {
        strcat(fixed_filename, "_fixed");
    }

    if (verbose) mexPrintf("Saving repaired header to: %s\n", fixed_filename);

    FILE *orig = fopen(edfFN, "rb");
    FILE *fixed = fopen(fixed_filename, "wb");
    if (!orig || !fixed) {
        if (verbose) mexPrintf("Error opening files for repair.\n");
        if (orig) fclose(orig);
        if (fixed) fclose(fixed);
        return;
    }

    fseek(orig, 0, SEEK_END);
    long size = ftell(orig);
    fseek(orig, 0, SEEK_SET);

    unsigned char *data = (unsigned char *)malloc(size);
    if (!data) { fclose(orig); fclose(fixed); return; }
    fread(data, 1, size, orig);

    char rec_str[16];
    sprintf(rec_str, "%-8d", num_records);
    for (int i = 0; i < 8; i++) data[236 + i] = (i < (int)strlen(rec_str)) ? rec_str[i] : ' ';

    fwrite(data, 1, size, fixed);
    free(data);
    fclose(orig);
    fclose(fixed);

    if (verbose) mexPrintf("Header repaired: num_data_records = %d\n", num_records);
}

/* --------------------------------------------------------------- */
/* Forward declarations                                            */
/* --------------------------------------------------------------- */
char *strtrim(char *str);

/* --------------------------------------------------------------- */
/* Annotation parsing                                              */
/* --------------------------------------------------------------- */
void parse_tals(unsigned char *data, int len, Annotation **annots, int *num_annots, int *capacity, int verbose) {
    int idx = 0;

    while (idx < len) {
        /* Skip leading zeros */
        while (idx < len && data[idx] == 0) idx++;
        if (idx >= len) break;

        /* Check for start of TAL (+ or -) */
        if (data[idx] != 43 && data[idx] != 45) {  /* + or - */
            idx++;
            continue;
        }

        int tal_start = idx;
        
        /* Read onset (from +/- up to first byte 20) */
        int onsetStart = idx;
        while (idx < len && data[idx] != 20) idx++;
        if (idx >= len) break;

        char onsetStr[256];
        int onsetLen = idx - onsetStart;
        if (onsetLen >= (int)sizeof(onsetStr)) onsetLen = sizeof(onsetStr) - 1;
        memcpy(onsetStr, &data[onsetStart], onsetLen);
        onsetStr[onsetLen] = '\0';
        
        double onset = atof(onsetStr);
        
        if (verbose) {
            mexPrintf("  TAL found at offset %d: onset_str='%s' -> %f\n", tal_start, onsetStr, onset);
        }

        idx++;  /* skip byte 20 */

        /* Read annotation texts (separated by byte 20, TAL ends with 0 or next +/-) */
        int text_count = 0;
        mxArray *texts_cell = mxCreateCellMatrix(1, 0);

        while (idx < len && data[idx] != 0 && data[idx] != 43 && data[idx] != 45) {
            /* Skip empty byte 20 separators */
            if (data[idx] == 20) {
                idx++;
                continue;
            }
            
            int textStart = idx;
            while (idx < len && data[idx] != 20 && data[idx] != 0 && data[idx] != 43 && data[idx] != 45) {
                idx++;
            }

            if (idx > textStart) {
                int textLen = idx - textStart;
                char text[512];
                if (textLen >= (int)sizeof(text)) textLen = sizeof(text) - 1;
                memcpy(text, &data[textStart], textLen);
                text[textLen] = '\0';

                if (verbose) {
                    mexPrintf("    Text %d: '%s'\n", text_count, text);
                }

                /* Expand cell array and add string */
                mxArray *new_cell = mxCreateCellMatrix(1, text_count + 1);
                for (int i = 0; i < text_count; i++) {
                    mxSetCell(new_cell, i, mxGetCell(texts_cell, i));
                }
                mxSetCell(new_cell, text_count, mxCreateString(text));
                mxDestroyArray(texts_cell);
                texts_cell = new_cell;
                text_count++;
            }
        }

        /* Only store TALs with actual text (skip time-keeping annotations) */
        if (text_count > 0) {
            if (*num_annots >= *capacity) {
                *capacity *= 2;
                Annotation *temp = (Annotation *)mxRealloc(*annots, (*capacity) * sizeof(Annotation));
                *annots = temp;
            }
            (*annots)[*num_annots].onset = onset;
            (*annots)[*num_annots].texts = texts_cell;
            (*num_annots)++;
            
            if (verbose) {
                mexPrintf("    Stored annotation %d with %d texts\n", *num_annots, text_count);
            }
        } else {
            mxDestroyArray(texts_cell);
            if (verbose) {
                mexPrintf("    Skipped (no text content)\n");
            }
        }
    }
}

void extract_annotations(const char *edfFN, EDF_Header *header, Signal_Header *sig_headers,
                         unsigned char *raw_bytes, long total_bytes_read, int verbose,
                         Annotation **annots, int *num_annots) {
    *annots = NULL;
    *num_annots = 0;
    int capacity = 100;  /* Initial capacity */

    /* Find annotation signal */
    int annotIdx = -1;
    for (int i = 0; i < header->num_signals; i++) {
        char trimmed[256];
        strcpy(trimmed, sig_headers[i].signal_labels);
        int si = 0, sj;
        while (trimmed[si] == ' ') si++;
        sj = 0;
        while (trimmed[si] != '\0') trimmed[sj++] = trimmed[si++];
        while (sj > 0 && trimmed[sj-1] == ' ') sj--;
        trimmed[sj] = '\0';
        
        if (verbose) {
            mexPrintf("Signal %d: '%s' (trimmed: '%s')\n", i, sig_headers[i].signal_labels, trimmed);
        }
        
        if (strcmp(trimmed, "EDF Annotations") == 0) {
            annotIdx = i;
            if (verbose) mexPrintf("Found EDF Annotations at signal index %d\n", annotIdx);
            break;
        }
    }

    if (annotIdx < 0) {
        if (verbose) mexPrintf("No EDF Annotations signal found. Num signals: %d\n", header->num_signals);
        return;  /* No annotation signal */
    }

    /* Calculate byte position and size of annotation signal per record */
    int annotBytePos = 0;
    for (int i = 0; i < annotIdx; i++) {
        annotBytePos += sig_headers[i].samples_in_record * 2;
    }

    int annotBytesPerRecord = sig_headers[annotIdx].samples_in_record * 2;
    int totalBytesPerRecord = 0;
    for (int i = 0; i < header->num_signals; i++) {
        totalBytesPerRecord += sig_headers[i].samples_in_record * 2;
    }

    if (verbose) {
        mexPrintf("Annotation signal info:\n");
        mexPrintf("  Byte position in record: %d\n", annotBytePos);
        mexPrintf("  Bytes per record: %d\n", annotBytesPerRecord);
        mexPrintf("  Total bytes per record: %d\n", totalBytesPerRecord);
        mexPrintf("  Samples in annot record: %d\n", sig_headers[annotIdx].samples_in_record);
    }

    /* Parse all records for annotations - but only scan first few and last */
    int num_records = total_bytes_read / totalBytesPerRecord;
    if (verbose) mexPrintf("Total records to scan: %d\n", num_records);
    
    unsigned char *data_ptr = raw_bytes;

    /* Allocate once at the beginning */
    *annots = (Annotation *)mxMalloc(capacity * sizeof(Annotation));
    *num_annots = 0;

    /* Scan all records but limit verbose output to first/last */
    int max_verbose_records = 5;
    for (int r = 0; r < num_records; r++) {
        unsigned char *record_ptr = data_ptr + r * totalBytesPerRecord;
        unsigned char *annotData = record_ptr + annotBytePos;

        if (verbose && (r < max_verbose_records || r >= (num_records - max_verbose_records))) {
            mexPrintf("Record %d annotation bytes (first 100): ", r);
            for (int i = 0; i < (annotBytesPerRecord < 100 ? annotBytesPerRecord : 100); i++) {
                if (annotData[i] >= 32 && annotData[i] < 127) {
                    mexPrintf("%c", annotData[i]);
                } else {
                    mexPrintf("[%d]", annotData[i]);
                }
            }
            mexPrintf("\n");
        } else if (verbose && r == max_verbose_records) {
            mexPrintf("... (skipping records %d to %d) ...\n", max_verbose_records, num_records - max_verbose_records - 1);
        }

        parse_tals(annotData, annotBytesPerRecord, annots, num_annots, &capacity, (r < 2 || r >= (num_records - 2)) ? verbose : 0);
    }
    
    if (verbose) mexPrintf("Total annotations parsed: %d\n", *num_annots);
}



/* --------------------------------------------------------------- */
/* Main mexFunction                                                */
/* --------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) mexErrMsgIdAndTxt("read_EDF_mex:InvalidInput", "Filename required.");

    /* ---- Parse inputs ---- */
    char edfFN[4096];
    if (!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("read_EDF_mex:InvalidInput", "First argument must be filename string.");
    mxGetString(prhs[0], edfFN, sizeof(edfFN));

    const mxArray *signal_labels_input = (nrhs >= 2) ? prhs[1] : NULL;
    const mxArray *epochs_input       = (nrhs >= 3) ? prhs[2] : NULL;
    int verbose = 1;
    if (nrhs >= 4) verbose = (mxGetScalar(prhs[3]) != 0.0);
    int repair_header = 0;
    if (nrhs >= 5) repair_header = (mxGetScalar(prhs[4]) != 0.0);

    /* ---- Open file ---- */
    FILE *fid = fopen(edfFN, "rb");
    if (!fid) mexErrMsgIdAndTxt("read_EDF_mex:FileError", "Cannot open file.");

    /* ---- Read headers ---- */
    EDF_Header header;
    if (!read_edf_header(fid, &header)) {
        fclose(fid);
        mexErrMsgIdAndTxt("read_EDF_mex:FileError", "Failed to read main header.");
    }

    Signal_Header *sig_headers = (Signal_Header *)mxCalloc(header.num_signals, sizeof(Signal_Header));
    if (!read_signal_headers(fid, sig_headers, header.num_signals)) {
        mxFree(sig_headers);
        fclose(fid);
        mexErrMsgIdAndTxt("read_EDF_mex:FileError", "Failed to read signal headers.");
    }

    /* Sampling frequencies */
    for (int i = 0; i < header.num_signals; i++) {
        sig_headers[i].sampling_frequency = sig_headers[i].samples_in_record / header.data_record_duration;
    }

    /* Total samples per record */
    int total_samples_per_record = 0;
    for (int i = 0; i < header.num_signals; i++) {
        total_samples_per_record += sig_headers[i].samples_in_record;
    }

    /* Fix invalid num_data_records */
    int num_records = header.num_data_records;
    if (num_records <= 0) {
        fseek(fid, 0, SEEK_END);
        long file_size = ftell(fid);
        long data_bytes = file_size - header.num_header_bytes;
        long bytes_per_record = (long)total_samples_per_record * 2;
        num_records = (int)(data_bytes / bytes_per_record);
        if (verbose) {
            mexWarnMsgIdAndTxt("read_EDF_mex:InvalidRecordCount",
                              "Invalid num_data_records (%d). Auto-detected: %d", header.num_data_records, num_records);
        }
    }

    if (verbose) {
        mexPrintf("=== EDF Summary ===\n");
        mexPrintf("Version: %s | Signals: %d | Records: %d | Duration: %.3f s\n",
                  header.edf_ver, header.num_signals, num_records, header.data_record_duration);
    }

    if (repair_header && header.num_data_records <= 0) {
        save_repaired_header(edfFN, num_records, verbose);
    }

    /* ---- Read raw data ---- */
    fseek(fid, 0, SEEK_END);
    long file_size = ftell(fid);
    fseek(fid, header.num_header_bytes, SEEK_SET);
    long total_bytes_expected = (long)total_samples_per_record * num_records * 2;
    unsigned char *raw_bytes = (unsigned char *)mxMalloc(total_bytes_expected);
    size_t bytes_read = fread(raw_bytes, 1, total_bytes_expected, fid);
    if (bytes_read != (size_t)total_bytes_expected) {
        num_records = (int)(bytes_read / (total_samples_per_record * 2));
        if (verbose) mexPrintf("Adjusted num_records to %d based on file size.\n", num_records);
    }
    fclose(fid);

    size_t total_samples = bytes_read / 2;
    int16_t *raw = (int16_t *)mxMalloc(total_samples * sizeof(int16_t));

    if (is_little_endian()) {
        memcpy(raw, raw_bytes, bytes_read);
    } else {
        for (size_t i = 0; i < total_samples; i++) {
            raw[i] = (int16_t)(raw_bytes[2*i] | (raw_bytes[2*i+1] << 8));
        }
    }

    /* ---- Determine signals to load ---- */
    int num_signals_to_load = header.num_signals;
    int *signal_indices = (int *)mxCalloc(header.num_signals, sizeof(int));

    if (signal_labels_input && mxIsCell(signal_labels_input) && mxGetNumberOfElements(signal_labels_input) > 0) {
        num_signals_to_load = 0;
        int nlabels = (int)mxGetNumberOfElements(signal_labels_input);
        for (int i = 0; i < nlabels; i++) {
            mxArray *cell = mxGetCell(signal_labels_input, i);
            if (!cell) continue;
            char label[256];
            mxGetString(cell, label, sizeof(label));
            int idx = signal_index_from_label(sig_headers, header.num_signals, label);
            if (idx >= 0) signal_indices[num_signals_to_load++] = idx;
        }
    } else {
        for (int i = 0; i < header.num_signals; i++) signal_indices[i] = i;
    }

    /* ---- Epoch selection ---- */
    int start_epoch = 0, end_epoch = num_records;
    if (epochs_input && mxGetNumberOfElements(epochs_input) >= 2) {
        double *e = mxGetPr(epochs_input);
        start_epoch = (int)e[0];
        end_epoch   = (int)e[1];
    }
    if (start_epoch < 0) start_epoch = 0;
    if (end_epoch > num_records) end_epoch = num_records;
    int num_epochs_to_load = end_epoch - start_epoch;

    if (verbose) {
        mexPrintf("Loading %d signals, epochs %d:%d (%d epochs)\n",
                  num_signals_to_load, start_epoch, end_epoch-1, num_epochs_to_load);
    }

    /* ---- Output 1: Main header ---- */
    if (nlhs >= 1) {
        const char *fields[] = {"edf_ver","patient_id","local_rec_id","recording_startdate",
                                "recording_starttime","num_header_bytes","num_data_records",
                                "data_record_duration","num_signals"};
        plhs[0] = mxCreateStructMatrix(1,1,9,fields);
        mxSetField(plhs[0],0,"edf_ver",               mxCreateString(header.edf_ver));
        mxSetField(plhs[0],0,"patient_id",            mxCreateString(header.patient_id));
        mxSetField(plhs[0],0,"local_rec_id",          mxCreateString(header.local_rec_id));
        mxSetField(plhs[0],0,"recording_startdate",   mxCreateString(header.recording_startdate));
        mxSetField(plhs[0],0,"recording_starttime",   mxCreateString(header.recording_starttime));
        mxSetField(plhs[0],0,"num_header_bytes",      mxCreateDoubleScalar(header.num_header_bytes));
        mxSetField(plhs[0],0,"num_data_records",      mxCreateDoubleScalar(num_records));
        mxSetField(plhs[0],0,"data_record_duration",  mxCreateDoubleScalar(header.data_record_duration));
        mxSetField(plhs[0],0,"num_signals",           mxCreateDoubleScalar(header.num_signals));
    }

    /* ---- Output 2: Signal headers ---- */
    if (nlhs >= 2) {
        const char *fields[] = {"signal_labels","transducer_type","physical_dimension",
                                "physical_min","physical_max","digital_min","digital_max",
                                "prefiltering","samples_in_record","sampling_frequency"};
        plhs[1] = mxCreateStructMatrix(1, num_signals_to_load, 10, fields);
        for (int i = 0; i < num_signals_to_load; i++) {
            int s = signal_indices[i];
            mxSetField(plhs[1],i,"signal_labels",       mxCreateString(sig_headers[s].signal_labels));
            mxSetField(plhs[1],i,"transducer_type",     mxCreateString(sig_headers[s].transducer_type));
            mxSetField(plhs[1],i,"physical_dimension",  mxCreateString(sig_headers[s].physical_dimension));
            mxSetField(plhs[1],i,"physical_min",        mxCreateDoubleScalar(sig_headers[s].physical_min));
            mxSetField(plhs[1],i,"physical_max",        mxCreateDoubleScalar(sig_headers[s].physical_max));
            mxSetField(plhs[1],i,"digital_min",         mxCreateDoubleScalar(sig_headers[s].digital_min));
            mxSetField(plhs[1],i,"digital_max",         mxCreateDoubleScalar(sig_headers[s].digital_max));
            mxSetField(plhs[1],i,"prefiltering",        mxCreateString(sig_headers[s].prefiltering));
            mxSetField(plhs[1],i,"samples_in_record",   mxCreateDoubleScalar(sig_headers[s].samples_in_record));
            mxSetField(plhs[1],i,"sampling_frequency",  mxCreateDoubleScalar(sig_headers[s].sampling_frequency));
        }
    }

    /* ---- Output 3: Signal data (serial conversion) ---- */
    if (nlhs >= 3) {
        plhs[2] = mxCreateCellMatrix(1, num_signals_to_load);

        for (int i = 0; i < num_signals_to_load; i++) {
            int s = signal_indices[i];
            int samples_per_record = sig_headers[s].samples_in_record;
            int ns = samples_per_record * num_epochs_to_load;
            
            mxArray *sig_array = mxCreateDoubleMatrix(1, ns, mxREAL);
            double *sig_data = mxGetPr(sig_array);
            mxSetCell(plhs[2], i, sig_array);

            double dig_range = sig_headers[s].digital_max - sig_headers[s].digital_min;
            double scale, offset;
            if (dig_range != 0.0) {
                scale  = (sig_headers[s].physical_max - sig_headers[s].physical_min) / dig_range;
                offset = sig_headers[s].physical_min - sig_headers[s].digital_min * scale;
            } else {
                scale  = 0.0;
                offset = sig_headers[s].physical_min;
            }

            /* Calculate base offset for this signal within a record */
            int base_offset = 0;
            for (int j = 0; j < s; j++) {
                base_offset += sig_headers[j].samples_in_record;
            }

            for (int r = 0; r < num_epochs_to_load; r++) {
                int rec_idx = start_epoch + r;
                long base_idx = (long)rec_idx * total_samples_per_record + base_offset;
                for (int samp = 0; samp < samples_per_record; samp++) {
                    long idx = base_idx + samp;
                    sig_data[r * samples_per_record + samp] = raw[idx] * scale + offset;
                }
            }
        }
    }

    /* ---- Output 4: Annotations ---- */
    if (nlhs >= 4) {
        Annotation *annots = NULL;
        int num_annots = 0;
        
        if (verbose) mexPrintf("\n=== Extracting Annotations ===\n");
        
        extract_annotations(edfFN, &header, sig_headers, raw_bytes, bytes_read + header.num_header_bytes, 
                           verbose, &annots, &num_annots);

        if (num_annots > 0) {
            const char *ann_fields[] = {"onset", "text"};
            plhs[3] = mxCreateStructMatrix(1, num_annots, 2, ann_fields);

            for (int i = 0; i < num_annots; i++) {
                mxSetField(plhs[3], i, "onset", mxCreateDoubleScalar(annots[i].onset));
                mxSetField(plhs[3], i, "text", annots[i].texts);
            }
            
            if (verbose) mexPrintf("Successfully created annotations struct with %d entries\n", num_annots);
            
            mxFree(annots);
        } else {
            if (verbose) mexPrintf("No annotations found, creating empty struct\n");
            plhs[3] = mxCreateStructMatrix(1, 0, 0, NULL);
        }
    }

    /* ---- Cleanup ---- */
    mxFree(raw);
    mxFree(raw_bytes);
    mxFree(signal_indices);
    mxFree(sig_headers);
}

/* Implementation of strtrim */
char *strtrim(char *str) {
    static char buffer[256];
    strcpy(buffer, str);
    int i = 0;
    while (buffer[i] == ' ') i++;
    int j = 0;
    while (buffer[i] != '\0') buffer[j++] = buffer[i++];
    while (j > 0 && buffer[j-1] == ' ') j--;
    buffer[j] = '\0';
    return buffer;
}