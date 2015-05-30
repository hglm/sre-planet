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

class ElevationScale {
public :
    // The spherical radius of the sea level or reference radius.
    float sea_level_radius;
    // Scale factor with which to multiply the 8-bit or 16-bit elevation value
    // normalized to [0, 1], to obtain the real proportional height.
    float scale_factor;
    // The elevation map value, normalized to [0, 1], corresponding to sea level.
    float normalized_elevation_map_sea_level;
};

extern const ElevationScale planet_elevation_scale[];

enum {
    MESH_FLAG_SPHERICAL_MODEL = 0x1,
    MESH_FLAG_NEW_CACHE = 0x2,
    MESH_FLAG_NO_CACHE = 0x4,
    MESH_FLAG_VERTEX_COLORS = 0x8
};

void CreatePlanetMesh(sreScene *scene,  const char *heightmap_filename,
    sreTexture *heightmap_texture, int detail_division_factor,
    sreTexture *color_map_texture, int color_map_horizontal_subdivisions,
    int color_map_vertical_subdivisions, unsigned int mesh_flags, int& nu_models,
    sreModel **&mesh_model);

void __attribute__ ((noreturn)) PlanetApplicationError(const char *format, ...);
