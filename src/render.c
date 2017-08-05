#include "render.h"
#include <webp/encode.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int render_system(Renderer *r, System *s, const char *filename) {
	memset(r->pixels,0,sizeof(r->pixels));

	for(int i = 0; i < s->N; i++) {
		int ix = s->parts[i].x[0]/s->L*IMG_W;
		int iy = s->parts[i].x[1]/s->L*IMG_H;

		for(int c = 0; c < 3; c++)
			r->pixels[3*(ix+IMG_W*iy)+c]=255;
	}

	uint8_t *output;
	int size = WebPEncodeRGB(r->pixels,IMG_W,IMG_H,IMG_W*3, 90, &output);
	if(size == 0)
		return 1;

	FILE *f = fopen(filename, "wb");
	if(f == 0)
		return 1;

	fwrite(output,size,1,f);

	fclose(f);
	free(output);
	return 0;
}
