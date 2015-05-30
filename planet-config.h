/*

Copyright (c) 2014 Harm Hanemaaijer <fgenfb@yahoo.com>

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

*/

// Number of seconds for a complete rotation of the earth.
#define DEFAULT_DAY_INTERVAL 1000.0f
#define YEAR_INTERVAL (DEFAULT_DAY_INTERVAL * 365.0f)

#define SPHERE
//#define CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
// The size of submeshes. Maximum MESH_WIDTH + 2 and MESH_HEIGHT + 2.
#define SUB_MESH_WIDTH 200
#define SUB_MESH_HEIGHT 200
// LONGITUDE and LATITUDE, in degrees, define the center of the part of the world that is shown
// when using plane. With a spherical model, it is the starting position.
// Negative longitude is west.
#define LONGITUDE 0
#define LATITUDE 0
// Exaggeration scaling factor of the terrain height relative the real proportions.
#define TERRAIN_HEIGHT_EXAGGERATION 50.0f
// Extra scaling factor applied to planet radius and terrain.
#define EXTRA_PLANET_SCALE_FACTOR 1.0f

#define SHOW_SPACECRAFT
