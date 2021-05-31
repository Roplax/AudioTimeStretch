#include <stdio.h>
#include <string.h>
#include <time.h>

// #include <stdexcept>
// #include "RunParameters.h"
// #include "WavFile.h"
// #include "SoundTouch.h"
// #include "BPMDetect.h"

#define DEBUG  1
#define FLOAT_SAMPLES  1
// floating point samples
typedef float  SAMPLETYPE;
// 16bit integer sample type
//    typedef short SAMPLETYPE;


// Processing chunk size (size chosen to be divisible by 2, 4, 6, 8, 10, 12, 14, 16 channels ...)
#define BUFF_SIZE           6720

//following saved for speech
#define SETTING_SEQUENCE_MS  40
#define SETTING_SEEKWINDOW_MS 15
#define SETTING_OVERLAP_MS  8

// for Music
#define SETTING_SEQUENCE_MS  82
#define SETTING_SEEKWINDOW_MS 28
#define SETTING_OVERLAP_MS  12

static int wsola_channels, wsola_sampleRate;
static double wsola_newTempofactor;

static int sequenceMs, seekWindowLength, seekLength,seekWindowMs, overlapLength;

// Setup the parameters for processing

/// Sets new tempo control value. Normal tempo = 1.0, smaller values
/// represent slower tempo, larger faster tempo.
    
static void setup(int sampleRate, int channels, double newTempofactor )
{
 
    int intskip;
    
    wsola_channels=channels;
    wsola_sampleRate=sampleRate;
    wsola_newTempofactor=newTempofactor;
    

    // Calculate new sequence duration
    calcSeqParameters();
    
    overlapLength=calculateOverlapLength(SETTING_OVERLAP_MS);

    // Calculate ideal skip length (according to tempo value) 
    nominalSkip = wsola_newTempofactor* (seekWindowLength - overlapLength);
    intskip = (int)(nominalSkip + 0.5);

    // Calculate how many samples are needed in the 'inputBuffer' to 
    // process another batch of samples
    //sampleReq = max(intskip + overlapLength, seekWindowLength) + seekLength / 2;
    sampleReq = max(intskip + overlapLength, seekWindowLength) + seekLength;

    if (DEBUG){
        // print processing information

    #ifdef INTEGER_SAMPLES
            fprintf(stderr, "Uses 16bit integer sample type in processing.\n\n");
    #else
        #ifndef FLOAT_SAMPLES
            #error "Sampletype not defined"
        #endif
            fprintf(stderr, "Uses 32bit floating point sample type in processing.\n\n");
    #endif
            // print processing information only if outFileName given i.e. some processing will happen
            fprintf(stderr, "Processing the file with the following parameters:\n");
            fprintf(stderr, "  channels = %i %%\n", wsola_channels);
            fprintf(stderr, "  sampleRate = %i \n", wsola_sampleRate);
            fprintf(stderr, "  tempofactor  = %i %%\n",wsola_newTempofactor);
            fflush(stderr);
    }
}

// Processes the sound
int process(SAMPLETYPE *in_buff, int insize, SAMPLETYPE *out_buff, int outsize)
    {
        int ovlSkip, offset;
        int temp;
    
        
    
        // If tempo differs from the normal ('SCALE'), scan for the best overlapping
        // position
        offset = seekBestOverlapPositionQuick(in_buff);

        // Mix the samples in the 'inputBuffer' at position of 'offset' with the 
        // samples in 'midBuffer' using sliding overlapping
        // ... first partially overlap with the end of the previous sequence
        // (that's in 'midBuffer')
        
        overlap(outputBuffer.ptrEnd((uint)overlapLength), inputBuffer.ptrBegin(), (uint)offset);
        
        overlapStereo(pOutput, pInput + 2 * ovlPos);
        outputBuffer.putSamples((uint)overlapLength);

        // ... then copy sequence samples from 'inputBuffer' to output:

        // length of sequence
        temp = (seekWindowLength - 2 * overlapLength);

        // crosscheck that we don't have buffer overflow...
        if ((int)inputBuffer.numSamples() < (offset + temp + overlapLength * 2))
        {
            continue;    // just in case, shouldn't really happen
        }

        outputBuffer.putSamples(inputBuffer.ptrBegin() + channels * (offset + overlapLength), (uint)temp);

        // Copies the end of the current sequence from 'inputBuffer' to 
        // 'midBuffer' for being mixed with the beginning of the next 
        // processing sequence and so on
        assert((offset + temp + overlapLength * 2) <= (int)inputBuffer.numSamples());
        
        memcpy(pMidBuffer, inputBuffer.ptrBegin() + channels * (offset + temp + overlapLength), 
            channels * sizeof(SAMPLETYPE) * overlapLength);

        // Remove the processed samples from the input buffer. Update
        // the difference between integer & nominal skip step to 'skipFract'
        // in order to prevent the error from accumulating over time.
        skipFract += nominalSkip;   // real skip size
        ovlSkip = (int)skipFract;   // rounded to integer skip
        skipFract -= ovlSkip;       // maintain the fraction part, i.e. real vs. integer skip
        inputBuffer.receiveSamples((uint)ovlSkip);
    
        return(0);
    }



// Quick seek algorithm for improved runtime-performance: First roughly scans through the 
// correlation area, and then scan surroundings of two best preliminary correlation candidates
// with improved precision
//
// Based on testing:
// - This algorithm gives on average 99% as good match as the full algorith
// - this quick seek algorithm finds the best match on ~90% of cases
// - on those 10% of cases when this algorithm doesn't find best match, 
//   it still finds on average ~90% match vs. the best possible match
int seekBestOverlapPositionQuick(const SAMPLETYPE *refPos,  SAMPLETYPE *pMidBuffer)
    {
        
    #define _MIN(a, b)   (((a) < (b)) ? (a) : (b))
    #define SCANSTEP    16
    #define SCANWIND    8

    int bestOffs;
    int i;
    int bestOffs2;
    float bestCorr, corr;
    float bestCorr2;
    double norm;

    // note: 'float' types used in this function in case that the platform would need to use software-fp

    bestCorr = FLT_MIN;
    bestOffs = SCANWIND;
    bestCorr2 = FLT_MIN;
    bestOffs2 = 0;

    int best = 0;

    // Scans for the best correlation value by testing each possible position
    // over the permitted range. Look for two best matches on the first pass to
    // increase possibility of ideal match.
    //
    // Begin from "SCANSTEP" instead of SCANWIND to make the calculation
    // catch the 'middlepoint' of seekLength vector as that's the a-priori 
    // expected best match position
    //
    // Roughly:
    // - 15% of cases find best result directly on the first round,
    // - 75% cases find better match on 2nd round around the best match from 1st round
    // - 10% cases find better match on 2nd round around the 2nd-best-match from 1st round
    for (i = SCANSTEP; i < seekLength - SCANWIND - 1; i += SCANSTEP)
    {
        // Calculates correlation value for the mixing position corresponding
        // to 'i'
        corr = (float)calcCrossCorr(refPos + wsola_channels*i, pMidBuffer, norm);
        // heuristic rule to slightly favour values close to mid of the seek range
        float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
        corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

        // Checks for the highest correlation value
        if (corr > bestCorr)
        {
            // found new best match. keep the previous best as 2nd best match
            bestCorr2 = bestCorr;
            bestOffs2 = bestOffs;
            bestCorr = corr;
            bestOffs = i;
        }
        else if (corr > bestCorr2)
        {
            // not new best, but still new 2nd best match
            bestCorr2 = corr;
            bestOffs2 = i;
        }
    }

    // Scans surroundings of the found best match with small stepping
    int end = _MIN(bestOffs + SCANWIND + 1, seekLength);
    for (i = bestOffs - SCANWIND; i < end; i++)
    {
        if (i == bestOffs) continue;    // this offset already calculated, thus skip

        // Calculates correlation value for the mixing position corresponding
        // to 'i'
        corr = (float)calcCrossCorr(refPos + wsola_channels*i, pMidBuffer, norm);
        // heuristic rule to slightly favour values close to mid of the range
        float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
        corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

        // Checks for the highest correlation value
        if (corr > bestCorr)
        {
            bestCorr = corr;
            bestOffs = i;
            best = 1;
        }
    }

    // Scans surroundings of the 2nd best match with small stepping
    end = _MIN(bestOffs2 + SCANWIND + 1, seekLength);
    for (i = bestOffs2 - SCANWIND; i < end; i++)
    {
        if (i == bestOffs2) continue;    // this offset already calculated, thus skip

        // Calculates correlation value for the mixing position corresponding
        // to 'i'
        corr = (float)calcCrossCorr(refPos + wsola_channels*i, pMidBuffer, norm);
        // heuristic rule to slightly favour values close to mid of the range
        float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
        corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

        // Checks for the highest correlation value
        if (corr > bestCorr)
        {
            bestCorr = corr;
            bestOffs = i;
            best = 2;
        }
    }

    // clear cross correlation routine state if necessary (is so e.g. in MMX routines).
    clearCrossCorrState();

    #ifdef INTEGER_SAMPLES
        adaptNormalizer();
    #endif

    return bestOffs;
    }

// Calculates processing sequence length according to tempo setting
void calcSeqParameters()
{
    // Adjust tempo param according to tempo, so that variating processing sequence length is used
    // at varius tempo settings, between the given low...top limits
    #define AUTOSEQ_TEMPO_LOW   0.5     // auto setting low tempo range (-50%)
    #define AUTOSEQ_TEMPO_TOP   2.0     // auto setting top tempo range (+100%)

    // sequence-ms setting values at above low & top tempo
    #define AUTOSEQ_AT_MIN      125.0
    #define AUTOSEQ_AT_MAX      50.0
    #define AUTOSEQ_K           ((AUTOSEQ_AT_MAX - AUTOSEQ_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW))
    #define AUTOSEQ_C           (AUTOSEQ_AT_MIN - (AUTOSEQ_K) * (AUTOSEQ_TEMPO_LOW))

    // seek-window-ms setting values at above low & top tempoq
    #define AUTOSEEK_AT_MIN     25.0
    #define AUTOSEEK_AT_MAX     15.0
    #define AUTOSEEK_K          ((AUTOSEEK_AT_MAX - AUTOSEEK_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW))
    #define AUTOSEEK_C          (AUTOSEEK_AT_MIN - (AUTOSEEK_K) * (AUTOSEQ_TEMPO_LOW))

    #define CHECK_LIMITS(x, mi, ma) (((x) < (mi)) ? (mi) : (((x) > (ma)) ? (ma) : (x)))

    double seq, seek;
    
 
    seq = AUTOSEQ_C + AUTOSEQ_K * wsola_newTempofactor;
    seq = CHECK_LIMITS(seq, AUTOSEQ_AT_MAX, AUTOSEQ_AT_MIN);
    sequenceMs = (int)(seq + 0.5);

    seek = AUTOSEEK_C + AUTOSEEK_K * wsola_newTempofactor;
    seek = CHECK_LIMITS(seek, AUTOSEEK_AT_MAX, AUTOSEEK_AT_MIN);
    seekWindowMs = (int)(seek + 0.5);
    

    // Update seek window lengths
    seekWindowLength = (wsola_sampleRate * sequenceMs) / 1000;
    if (seekWindowLength < 2 * overlapLength) 
    {
        seekWindowLength = 2 * overlapLength;
    }
    seekLength = (wsola_sampleRate * seekWindowMs) / 1000;
}



//////////////////////////////////////////////////////////////////////////////
//
// Floating point arithmetics specific algorithm implementations.
//

#ifdef FLOAT_SAMPLES

// Overlaps samples in 'midBuffer' with the samples in 'pInput'
void overlapStereo(float *pOutput, const float *pInput) 
{
    int i;
    float fScale;
    float f1;
    float f2;

    fScale = 1.0f / (float)overlapLength;

    f1 = 0;
    f2 = 1.0f;

    for (i = 0; i < 2 * (int)overlapLength ; i += 2) 
    {
        pOutput[i + 0] = pInput[i + 0] * f1 + pMidBuffer[i + 0] * f2;
        pOutput[i + 1] = pInput[i + 1] * f1 + pMidBuffer[i + 1] * f2;

        f1 += fScale;
        f2 -= fScale;
    }
}


// Overlaps samples in 'midBuffer' with the samples in 'input'. 
void overlapMulti(float *pOutput, const float *pInput) 
{
    int i;
    float fScale;
    float f1;
    float f2;

    fScale = 1.0f / (float)overlapLength;

    f1 = 0;
    f2 = 1.0f;

    i=0;
    for (int i2 = 0; i2 < overlapLength; i2 ++)
    {
        // note: Could optimize this slightly by taking into account that always channels > 2
        for (int c = 0; c < wsola_channels; c ++)
        {
            pOutput[i] = pInput[i] * f1 + pMidBuffer[i] * f2;
            i++;
        }
        f1 += fScale;
        f2 -= fScale;
    }
}


/// Calculates overlapInMsec period length in samples.
int calculateOverlapLength(int overlapInMsec)
{
    int newOvl;

    assert(overlapInMsec >= 0);
    newOvl = (wsola_sampleRate * overlapInMsec) / 1000;
    if (newOvl < 16) newOvl = 16;

    // must be divisible by 8
    newOvl -= newOvl % 8;

    return(newOvl);
}


/// Calculate cross-correlation
double calcCrossCorr(const float *mixingPos, const float *compare, double *anorm)
{
    double corr;
    double norm;
    int i;

    corr = norm = 0;
    // Same routine for stereo and mono. For Stereo, unroll by factor of 2.
    // For mono it's same routine yet unrollsd by factor of 4.
    for (i = 0; i < wsola_channels * overlapLength; i += 4) 
    {
        corr += mixingPos[i] * compare[i] +
                mixingPos[i + 1] * compare[i + 1];

        norm += mixingPos[i] * mixingPos[i] + 
                mixingPos[i + 1] * mixingPos[i + 1];

        // unroll the loop for better CPU efficiency:
        corr += mixingPos[i + 2] * compare[i + 2] +
                mixingPos[i + 3] * compare[i + 3];

        norm += mixingPos[i + 2] * mixingPos[i + 2] +
                mixingPos[i + 3] * mixingPos[i + 3];
    }

    anorm = norm;
    return corr / sqrt((norm < 1e-9 ? 1.0 : norm));
}


/// Update cross-correlation by accumulating "norm" coefficient by previously calculated value
double calcCrossCorrAccumulate(const float *mixingPos, const float *compare, double *norm)
{
    double corr;
    int i;

    corr = 0;

    // cancel first normalizer tap from previous round
    for (i = 1; i <= wsola_channels; i ++)
    {
        norm -= mixingPos[-i] * mixingPos[-i];
    }

    // Same routine for stereo and mono. For Stereo, unroll by factor of 2.
    // For mono it's same routine yet unrollsd by factor of 4.
    for (i = 0; i < wsola_channels * overlapLength; i += 4) 
    {
        corr += mixingPos[i] * compare[i] +
                mixingPos[i + 1] * compare[i + 1] +
                mixingPos[i + 2] * compare[i + 2] +
                mixingPos[i + 3] * compare[i + 3];
    }

    // update normalizer with last samples of this round
    for (int j = 0; j < wsola_channels; j ++)
    {
        i --;
        norm += mixingPos[i] * mixingPos[i];
    }

    return corr / sqrt((norm < 1e-9 ? 1.0 : norm));
}


#endif // FLOAT_SAMPLES


//////////////////////////////////////////////////////////////////////////////
//
// Integer arithmetics specific algorithm implementations.
//
//////////////////////////////////////////////////////////////////////////////

#ifdef INTEGER_SAMPLES

/// For integer algorithm: adapt normalization factor divider with music so that 
/// it'll not be pessimistically restrictive that can degrade quality on quieter sections
/// yet won't cause integer overflows either

void adaptNormalizer()
{
    // Do not adapt normalizer over too silent sequences to avoid averaging filter depleting to
    // too low values during pauses in music
    if ((maxnorm > 1000) || (maxnormf > 40000000))
    { 
        //norm averaging filter
        maxnormf = 0.9f * maxnormf + 0.1f * (float)maxnorm;

        if ((maxnorm > 800000000) && (overlapDividerBitsNorm < 16))
        {
            // large values, so increase divider
            overlapDividerBitsNorm++;
            if (maxnorm > 1600000000) overlapDividerBitsNorm++; // extra large value => extra increase
        }
        else if ((maxnormf < 1000000) && (overlapDividerBitsNorm > 0))
        {
            // extra small values, decrease divider
            overlapDividerBitsNorm--;
        }
    }

    maxnorm = 0;
}

// Overlaps samples in 'midBuffer' with the samples in 'input'. The 'Stereo' 
// version of the routine.
void overlapStereo(short *poutput, const short *input) const
{
    int i;
    short temp;
    int cnt2;

    for (i = 0; i < overlapLength ; i ++) 
    {
        temp = (short)(overlapLength - i);
        cnt2 = 2 * i;
        poutput[cnt2] = (input[cnt2] * i + pMidBuffer[cnt2] * temp )  / overlapLength;
        poutput[cnt2 + 1] = (input[cnt2 + 1] * i + pMidBuffer[cnt2 + 1] * temp ) / overlapLength;
    }
}


// Overlaps samples in 'midBuffer' with the samples in 'input'. The 'Multi'
// version of the routine.
void overlapMulti(SAMPLETYPE *poutput, const SAMPLETYPE *input) const
{
    SAMPLETYPE m1=(SAMPLETYPE)0;
    SAMPLETYPE m2;
    int i=0;

    for (m2 = (SAMPLETYPE)overlapLength; m2; m2 --)
    {
        for (int c = 0; c < wsola_channels; c ++)
        {
            poutput[i] = (input[i] * m1 + pMidBuffer[i] * m2)  / overlapLength;
            i++;
        }

        m1++;
    }
}

// Calculates the x having the closest 2^x value for the given value
static int _getClosest2Power(double value)
{
    return (int)(log(value) / log(2.0) + 0.5);
}


/// Calculates overlap period length in samples.
/// Integer version rounds overlap length to closest power of 2
/// for a divide scaling operation.
int calculateOverlapLength(int aoverlapMs)
{
    int newOvl;

    assert(aoverlapMs >= 0);

    // calculate overlap length so that it's power of 2 - thus it's easy to do
    // integer division by right-shifting. Term "-1" at end is to account for 
    // the extra most significatnt bit left unused in result by signed multiplication 
    overlapDividerBitsPure = _getClosest2Power((wsola_sampleRate * aoverlapMs) / 1000.0) - 1;
    if (overlapDividerBitsPure > 9) overlapDividerBitsPure = 9;
    if (overlapDividerBitsPure < 3) overlapDividerBitsPure = 3;
    newOvl = (int)pow(2.0, (int)overlapDividerBitsPure + 1);    // +1 => account for -1 above

    

    //overlapDividerBitsNorm = overlapDividerBitsPure;

    // calculate sloping divider so that crosscorrelation operation won't 
    // overflow 32-bit register. Max. sum of the crosscorrelation sum without 
    // divider would be 2^30*(N^3-N)/3, where N = overlap length
    //slopingDivider = (newOvl * newOvl - 1) / 3;
    
    return(newOvl);
}


double calcCrossCorr(const short *mixingPos, const short *compare, double &norm)
{
    long corr;
    unsigned long lnorm;
    int i;

    corr = lnorm = 0;
    // Same routine for stereo and mono. For stereo, unroll loop for better
    // efficiency and gives slightly better resolution against rounding. 
    // For mono it same routine, just  unrolls loop by factor of 4
    for (i = 0; i < wsola_channels * overlapLength; i += 4) 
    {
        corr += (mixingPos[i] * compare[i] + 
                 mixingPos[i + 1] * compare[i + 1]) >> overlapDividerBitsNorm;  // notice: do intermediate division here to avoid integer overflow
        corr += (mixingPos[i + 2] * compare[i + 2] + 
                mixingPos[i + 3] * compare[i + 3]) >> overlapDividerBitsNorm;
        lnorm += (mixingPos[i] * mixingPos[i] + 
                mixingPos[i + 1] * mixingPos[i + 1]) >> overlapDividerBitsNorm; // notice: do intermediate division here to avoid integer overflow
        lnorm += (mixingPos[i + 2] * mixingPos[i + 2] + 
                mixingPos[i + 3] * mixingPos[i + 3]) >> overlapDividerBitsNorm;
    }

    if (lnorm > maxnorm)
    {
        maxnorm = lnorm;
    }
    // Normalize result by dividing by sqrt(norm) - this step is easiest 
    // done using floating point operation
    norm = (double)lnorm;
    return (double)corr / sqrt((norm < 1e-9) ? 1.0 : norm);
}


/// Update cross-correlation by accumulating "norm" coefficient by previously calculated value
double calcCrossCorrAccumulate(const short *mixingPos, const short *compare, double &norm)
{
    long corr;
    unsigned long lnorm;
    int i;

    // cancel first normalizer tap from previous round
    lnorm = 0;
    for (i = 1; i <= wsola_channels; i ++)
    {
        lnorm -= (mixingPos[-i] * mixingPos[-i]) >> overlapDividerBitsNorm;
    }

    corr = 0;
    // Same routine for stereo and mono. For stereo, unroll loop for better
    // efficiency and gives slightly better resolution against rounding. 
    // For mono it same routine, just  unrolls loop by factor of 4
    for (i = 0; i < wsola_channels * overlapLength; i += 4) 
    {
        corr += (mixingPos[i] * compare[i] + 
                 mixingPos[i + 1] * compare[i + 1]) >> overlapDividerBitsNorm;  // notice: do intermediate division here to avoid integer overflow
        corr += (mixingPos[i + 2] * compare[i + 2] + 
                 mixingPos[i + 3] * compare[i + 3]) >> overlapDividerBitsNorm;
    }

    // update normalizer with last samples of this round
    for (int j = 0; j < wsola_channels; j ++)
    {
        i --;
        lnorm += (mixingPos[i] * mixingPos[i]) >> overlapDividerBitsNorm;
    }

    norm += (double)lnorm;
    if (norm > maxnorm)
    {
        maxnorm = (unsigned long)norm;
    }

    // Normalize result by dividing by sqrt(norm) - this step is easiest 
    // done using floating point operation
    return (double)corr / sqrt((norm < 1e-9) ? 1.0 : norm);
}

#endif // INTEGER_SAMPLES

