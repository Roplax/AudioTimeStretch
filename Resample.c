#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "wav.h"

int main(int argc, char *argv[]) {
	char* filename = argv[1];
	int stretchFactor = atoi(argv[2]);
	if (stretchFactor < 1) {
		printf("Invalid stretch factor\n");
		return 1;
	}
	struct wav_info info;

	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Error reading file\n");
		return 1;
	}

	char* outputName;
	asprintf(&outputName, "STRETCH_%s", filename);
	FILE* op = fopen(outputName, "wb");
	if (!op) {
		printf("Error writing file\n");
		return 1;
	}


	read_wav_info(&info, fp);

	int originalSampleNum = info.num_samples;
	info.num_samples *= stretchFactor;
	write_wav_hdr(&info, op);

	int_fast32_t sample[info.num_channels];
	for (int i = 0; i < originalSampleNum; i++) {
		for (int j = 0; j < info.num_channels; j++) {
			fread(&sample[j], info.bits_per_sample/8, 1, fp);
		}
		for (int k = 0; k < stretchFactor; k++) {
			write_sample(&info, op, sample);
		}
	}
	fclose(fp);
	fclose(op);
	printf("Successfully stretched file by given factor\n");
}
