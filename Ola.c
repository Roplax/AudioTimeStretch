#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "wav.h"

#define OVERLAP 0.5
#define WINDOWLENGTH 4096
#define MaxChannels 2
#define SIZECAP 4000
#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

//

//works better with making files shorter than it does with making them longer
//use stretchFactor < 1 for shorter files
//hopefully wsola will be better than ola in this regard


int main(int argc, char *argv[]) {
	char* filename = argv[1];
	float stretchFactor = atof(argv[2]);

	struct wav_info info;

	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Error Opening file\n");
		return 1;
	}

	char* outputName;
	asprintf(&outputName, "OLA_%s", filename);
	FILE* op = fopen(outputName, "wb");
	if (!op) {
		printf("Error writing file\n");
		return 1;
	}


	read_wav_info(&info, fp);

	int sizeInBits = info.num_samples * info.num_channels * info.bits_per_sample;

	int originalSampleNum = info.num_samples;

	int numChunks=originalSampleNum/SIZECAP; //should this be checked against max int?

	if (sizeInBits > SIZECAP) {
		printf("File is too big, using chunking to process.  Actual Data Size:%d \n",sizeInBits/8 );
		//return 1;
	} else if (stretchFactor > 3.99) {
		printf("stretchFactor is too big. Aborting \n");
		return 1;
	}


	info.num_samples *= stretchFactor;
	write_wav_hdr(&info, op);

	//initialize all the variables for the ola algorithm
	int windowOffset = (1 - OVERLAP) * WINDOWLENGTH;
	int gammaLength = (int) info.num_samples / windowOffset;
	int sigmaInterval = windowOffset / stretchFactor;
	int gamma[gammaLength];
	int sigma[gammaLength];
	float window[sigmaInterval];

	//populate gamma and sigma
	gamma[0] = sigma[0] = 0;
	for (int i = 1; i < gammaLength; i++) {
		gamma[i] = gamma[i-1] + windowOffset;
		sigma[i] = gamma[i] + sigmaInterval;
	}

	//populate window
	for (int i = 0; i < sigmaInterval; i++) {
		window[i] = pow(sin(M_PI * i / (sigmaInterval - 1)), 2);
	}

	//initialize input and output buffers
	printf("Initializing input and output buffers for %d samples and %d channels\n",originalSampleNum, info.num_channels);
	printf("Initializing allSamples\n");
	int_fast32_t allSamples[SIZECAP][MaxChannels];
	printf("Initializing outputSamples\n");
	int_fast32_t outputSamples[SIZECAP*8][MaxChannels];

        int SamplesTobeRead,SamplesToWrite, TotalSamplesRead;
        SamplesTobeRead=min(SIZECAP,originalSampleNum);

        //Read input and process in chunks
        for (int ichunk = 0; ichunk < numChunks+1; ichunk++) {

                SamplesToWrite=SamplesTobeRead*stretchFactor;

		//zero out outputSamples
		printf("Zeroing out outputSamples\n");        
		for (int i = 0; i < SamplesToWrite; i++) {
			for (int j = 0; j < info.num_channels; j++) {
				outputSamples[i][j] = 0;
			}
		}

		//read all samples into memory
		printf("Read %d samples into memory Of Total: %d \n",SamplesTobeRead, TotalSamplesRead );
		for (int i = 0; i < SamplesTobeRead; i++) {
			for (int j = 0; j < info.num_channels; j++) {
				fread(&allSamples[i][j], info.bits_per_sample/8, 1, fp);
			}
		}
		TotalSamplesRead += SamplesTobeRead;

		printf("Processing to produce outputSamples\n");  
		//read each frame, window it, and copy it to the output buffer
		for (int i = 0; i < gammaLength - 2; i++) { //loop through each element of sigma
			int_fast32_t frame[sigmaInterval][info.num_channels];

			for (int j = 0; j < sigmaInterval; j++) { //loop through each sample of the frame
				for (int k = 0; k < info.num_channels; k++) { //loop through each channel of the sample
					frame[j][k] = allSamples[sigma[i] + j][k] * window[j];
				}
			} 
			for (int l = 0; l < sigmaInterval; l++) { //loop through each sample of frame again
				for (int m = 0; m < info.num_channels; m++) { //to add each sample to output buffer
					outputSamples[gamma[i] + l][m] += frame[l][m]; //Segfault here where OutputSamples 4*SIZECAP 
				}
			}

		}
		printf("Writing %d outputSamples\n", SamplesToWrite);  
		//write all samples
		int_fast32_t sample[info.num_channels];
		for (int i = 0; i <SamplesToWrite; i++) {
			for (int j = 0; j < info.num_channels; j++) {
				sample[j] = outputSamples[i][j];
			}
			write_sample(&info, op, sample);
		}
		SamplesTobeRead=min(SIZECAP,(originalSampleNum-TotalSamplesRead));

        }


	fclose(fp);
	fclose(op);
	printf("Successfully stretched file by given factor\n");
}
