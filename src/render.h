#pragma once
#include <stdint.h>
#include "mcgas.h"
enum {
	IMG_W = 400,
	IMG_H = 400
};

typedef struct {
	uint8_t pixels[3*IMG_W*IMG_H];
} Renderer;

int render_system(Renderer *r, System *s, const char *filename);
