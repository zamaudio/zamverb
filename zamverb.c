#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define ZAMVERB_URI "http://zamaudio.com/lv2/zamverb"
#define MAXALLPS 4
#define MAXCOMBS 8
#define MAXBUFFER 2000
#define SPREAD 23
#define PI 3.1415926
#define BANDPASSF0 2000
#define BANDPASSQ 0.9

typedef enum {
	ZAMVERB_INPUTL  = 0,
	ZAMVERB_INPUTR  = 1,
	ZAMVERB_OUTPUTL = 2,
	ZAMVERB_OUTPUTR = 3,
	
	ZAMVERB_ROOMSIZE = 4,
	ZAMVERB_WET = 5
} PortIndex;

typedef struct {
	float* inputl;
	float* inputr;
	float* outputl;
	float* outputr;
	float* roomsize;
	float* wet;

	float allpassfeedback;
	float combfeedback;
	float combdamp1;
	float combdamp2;
	float filterbufL[MAXCOMBS];
	float filterbufR[MAXCOMBS];

	float bufcombL[MAXCOMBS][MAXBUFFER];
	float bufcombR[MAXCOMBS][MAXBUFFER];
	int combidxL[MAXCOMBS];
	int combidxR[MAXCOMBS];

	float bufallpassL[MAXALLPS][MAXBUFFER];
	float bufallpassR[MAXALLPS][MAXBUFFER];
	int allpassidxL[MAXALLPS];
	int allpassidxR[MAXALLPS];

	int combtuneL[MAXCOMBS];
	int combtuneR[MAXCOMBS];

	int allpasstuneL[MAXALLPS];
	int allpasstuneR[MAXALLPS];

	float lpx1L;
	float lpx2L;
	float lpy1L;
	float lpy2L;
	float lpx1R;
	float lpx2R;
	float lpy1R;
	float lpy2R;
	float lpa0;
	float lpa1;
	float lpa2;
	float lpb0;
	float lpb1;
	float lpb2;

	float srate;
} ZamVERB;



static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               bundle_path,
            const LV2_Feature* const* features)
{
	int i;
	ZamVERB* zamverb = (ZamVERB*)malloc(sizeof(ZamVERB));
	zamverb->srate = rate;

	zamverb->combtuneL[0] = 1116;
	zamverb->combtuneL[1] = 1188;
	zamverb->combtuneL[2] = 1277;
	zamverb->combtuneL[3] = 1356;
	zamverb->combtuneL[4] = 1422;
	zamverb->combtuneL[5] = 1491;
	zamverb->combtuneL[6] = 1557;
	zamverb->combtuneL[7] = 1617;
	zamverb->combtuneR[0] = 1116+SPREAD;
	zamverb->combtuneR[1] = 1188+SPREAD;
	zamverb->combtuneR[2] = 1277+SPREAD;
	zamverb->combtuneR[3] = 1356+SPREAD;
	zamverb->combtuneR[4] = 1422+SPREAD;
	zamverb->combtuneR[5] = 1491+SPREAD;
	zamverb->combtuneR[6] = 1557+SPREAD;
	zamverb->combtuneR[7] = 1617+SPREAD;

	zamverb->allpasstuneL[0] = 556;
	zamverb->allpasstuneL[1] = 441;
	zamverb->allpasstuneL[2] = 341;
	zamverb->allpasstuneL[3] = 225;
	zamverb->allpasstuneR[0] = 556+SPREAD;
	zamverb->allpasstuneR[1] = 441+SPREAD;
	zamverb->allpasstuneR[2] = 341+SPREAD;
	zamverb->allpasstuneR[3] = 224+SPREAD;

        zamverb->lpx1L = 0.f;
        zamverb->lpx2L = 0.f;
        zamverb->lpy1L = 0.f;
        zamverb->lpy2L = 0.f;
        zamverb->lpx1R = 0.f;
        zamverb->lpx2R = 0.f;
        zamverb->lpy1R = 0.f;
        zamverb->lpy2R = 0.f;

	for (i = 0; i < MAXCOMBS; i++) {
		zamverb->combidxL[i] = 0;
		zamverb->combidxR[i] = 0;
		zamverb->filterbufL[i] = 0.f;
		zamverb->filterbufR[i] = 0.f;
		zamverb->combtuneL[i] = (int)(zamverb->combtuneL[i] * rate / 44100);
		zamverb->combtuneR[i] = (int)(zamverb->combtuneR[i] * rate / 44100);
		memset(&zamverb->bufcombL[i], 0x0, MAXBUFFER*sizeof(float));
		memset(&zamverb->bufcombR[i], 0x0, MAXBUFFER*sizeof(float));
	}
	
        for (i = 0; i < MAXALLPS; i++) {
		zamverb->allpassidxL[i] = 0;
		zamverb->allpassidxR[i] = 0;
		zamverb->allpasstuneL[i] = (int)(zamverb->allpasstuneL[i] * rate / 44100);
		zamverb->allpasstuneR[i] = (int)(zamverb->allpasstuneR[i] * rate / 44100);
		memset(&zamverb->bufallpassL[i], 0x0, MAXBUFFER*sizeof(float));
		memset(&zamverb->bufallpassR[i], 0x0, MAXBUFFER*sizeof(float));
        }

	return (LV2_Handle)zamverb;
}

static void
connect_port(LV2_Handle instance,
             uint32_t   port,
             void*      data)
{
	ZamVERB* zamverb = (ZamVERB*)instance;

	switch ((PortIndex)port) {
	case ZAMVERB_INPUTL:
		zamverb->inputl = (float*)data;
		break;
	case ZAMVERB_INPUTR:
		zamverb->inputr = (float*)data;
		break;
	case ZAMVERB_OUTPUTL:
		zamverb->outputl = (float*)data;
		break;
	case ZAMVERB_OUTPUTR:
		zamverb->outputr = (float*)data;
		break;
	case ZAMVERB_ROOMSIZE:
		zamverb->roomsize = (float*)data;
		break;
	case ZAMVERB_WET:
		zamverb->wet = (float*)data;
		break;
	}
}

// Force already-denormal float value to zero
static inline void
sanitize_denormal(float *value) {
    if (!isnormal(*value)) {
        *value = 0.f;
    }
}

static inline int 
sign(float x) {
        return (x >= 0.f ? 1 : -1);
}

static inline float 
from_dB(float gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline float
to_dB(float g) {
        return (20.f*log10(g));
}

static void
activate(LV2_Handle instance)
{
}

static float
bandpass_procL(LV2_Handle instance, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float output;
	output = input*zamverb->lpb0 
		+ zamverb->lpx1L*zamverb->lpb1
		+ zamverb->lpx2L*zamverb->lpb2
		- zamverb->lpy1L*zamverb->lpa1
		- zamverb->lpy2L*zamverb->lpa2;

	sanitize_denormal(&output);

	zamverb->lpx2L = zamverb->lpx1L;
	zamverb->lpx1L = input;
	zamverb->lpy2L = zamverb->lpy1L;
	zamverb->lpy1L = output;

	return output;
}

static float
bandpass_procR(LV2_Handle instance, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float output;
	output = input*zamverb->lpb0 
		+ zamverb->lpx1R*zamverb->lpb1 
		+ zamverb->lpx2R*zamverb->lpb2
		- zamverb->lpy1R*zamverb->lpa1
		- zamverb->lpy2R*zamverb->lpa2;

	sanitize_denormal(&output);

	zamverb->lpx2R = zamverb->lpx1R;
	zamverb->lpx1R = input;
	zamverb->lpy2R = zamverb->lpy1R;
	zamverb->lpy1R = output;

	return output;
}

static float
allpass_procL(LV2_Handle instance, int fx, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float output;
	float bufout;
	int bx;
	float feedback;
	
	feedback = zamverb->allpassfeedback;
	bx = zamverb->allpassidxL[fx];
	bufout = (float) zamverb->bufallpassL[fx][bx];

	sanitize_denormal(&bufout);

	output = -input + bufout;
	zamverb->bufallpassL[fx][bx] = input + bufout*feedback;

	if (++zamverb->allpassidxL[fx] >= zamverb->allpasstuneL[fx])
		zamverb->allpassidxL[fx] = 0;

	return output;
}

static float
allpass_procR(LV2_Handle instance, int fx, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float output;
	float bufout;
	int bx;
	float feedback = zamverb->allpassfeedback;

	bx = zamverb->allpassidxR[fx];
	bufout = (float) zamverb->bufallpassR[fx][bx];

	sanitize_denormal(&bufout);

	output = -input + bufout;
	zamverb->bufallpassR[fx][bx] = input + bufout*feedback;

	if (++zamverb->allpassidxR[fx] >= zamverb->allpasstuneR[fx])
		zamverb->allpassidxR[fx] = 0;

	return output;
}

static float
comb_procL(LV2_Handle instance, int fx, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float damp1 = zamverb->combdamp1;
	float damp2 = zamverb->combdamp2;
	float feedback = zamverb->combfeedback;
	float output;
	int bx;

	bx = zamverb->combidxL[fx];
	output = zamverb->bufcombL[fx][bx];

	sanitize_denormal(&output);

	zamverb->filterbufL[fx] = output*damp2
		+ zamverb->filterbufL[fx]*damp1;
	
	sanitize_denormal(&zamverb->filterbufL[fx]);

	zamverb->bufcombL[fx][bx] = input 
		+ zamverb->filterbufL[fx]*feedback;
	
	if (++zamverb->combidxL[fx] >= zamverb->combtuneL[fx])
		zamverb->combidxL[fx] = 0;
	
	return output;
}

static float
comb_procR(LV2_Handle instance, int fx, float input)
{
	ZamVERB* zamverb = (ZamVERB*)instance;
	float damp1 = zamverb->combdamp1;
	float damp2 = zamverb->combdamp2;
	float feedback = zamverb->combfeedback;
	float output;
	int bx;

	bx = zamverb->combidxR[fx];
	output = zamverb->bufcombR[fx][bx];

	sanitize_denormal(&output);

	zamverb->filterbufR[fx] = output*damp2
		+ zamverb->filterbufR[fx]*damp1;
	
	sanitize_denormal(&zamverb->filterbufR[fx]);

	zamverb->bufcombR[fx][bx] = input 
		+ zamverb->filterbufR[fx]*feedback;
	
	if (++zamverb->combidxR[fx] >= zamverb->combtuneR[fx])
		zamverb->combidxR[fx] = 0;
	
	return output;
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	ZamVERB* zamverb = (ZamVERB*)instance;

	const float* const inputl  = zamverb->inputl;
	const float* const inputr  = zamverb->inputr;
	float* const       outputl = zamverb->outputl;
	float* const       outputr = zamverb->outputr;

	float roomsize = (*(zamverb->roomsize)-0.7)/0.28;
	float damp = 0.05;
	float wet = *(zamverb->wet)/2.f;
	float dry = (1.f - *(zamverb->wet)/2.f);
	
	zamverb->combdamp1 = damp;
	zamverb->combdamp2 = 1.f - damp;
	zamverb->combfeedback = roomsize;
	zamverb->allpassfeedback = roomsize;

	float omega = 2.0 * PI * BANDPASSF0 / zamverb->srate;
	float alpha = sin(omega) / (2.0 * BANDPASSQ);

	zamverb->lpa0 = 1.0 + alpha;
	zamverb->lpa1 = -2.0*cos(omega) / zamverb->lpa0;
	zamverb->lpa2 = (1.0 - alpha) / zamverb->lpa0;

	zamverb->lpb0 = alpha / zamverb->lpa0;
	zamverb->lpb1 = 0.0;
	zamverb->lpb2 = -alpha / zamverb->lpa0;
	
	uint32_t pos;
	for (pos = 0; pos < n_samples; pos++) {
		float inl = inputl[pos];
		float inr = inputr[pos];
		float outl = 0.f;
		float outr = 0.f;
		int i;

		// Comb filters
		for (i = 0; i < MAXCOMBS; ++i) {
			outl += comb_procL(zamverb, i, inl);
			outr += comb_procR(zamverb, i, inr);
		}

		// Allpass filters
		for (i = 0; i < MAXALLPS; ++i) {
			outl = allpass_procL(zamverb, i, outl);
			outr = allpass_procR(zamverb, i, outr);
		}

		// Bandpass filter
		outl = bandpass_procL(zamverb, outl);
		outr = bandpass_procR(zamverb, outr);

		// Mix wet/dry
		outputl[pos] = (outl*wet/8.0 + inl*dry)/2.f;
		outputr[pos] = (outr*wet/8.0 + inr*dry)/2.f;
	}
}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
	free(instance);
}

const void*
extension_data(const char* uri)
{
	return NULL;
}

static const LV2_Descriptor descriptor = {
	ZAMVERB_URI,
	instantiate,
	connect_port,
	activate,
	run,
	deactivate,
	cleanup,
	extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
	switch (index) {
	case 0:
		return &descriptor;
	default:
		return NULL;
	}
}
