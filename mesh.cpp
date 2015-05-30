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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "sre.h"
#include "sreBackend.h"

#include "planet-config.h"
#include "sre-planet.h"

#define MESH_INDEX(x, y) (y * WIDTH_IN_MESH_MODELS + x)

static const char *planet_elevation_map_filename;
static sreTexture *planet_elevation_map_texture;
static int planet_elevation_map_detail_scale_down;
// Planet terrain color texture variables used when MESH_FLAG_VERTEX_COLORS is defined.
static sreTexture *planet_color_map_texture;
static int color_map_horizontal_subdivisions, color_map_vertical_subdivisions;
static unsigned int planet_mesh_flags;
static sreModel **planet_mesh_model;
static int MESH_WIDTH, MESH_HEIGHT;
static int MESH_SAMPLE_AREA_WIDTH, MESH_SAMPLE_AREA_HEIGHT;
static int WIDTH_IN_MESH_MODELS, HEIGHT_IN_MESH_MODELS;
static int X_OFFSET, Y_OFFSET;
// Zoom must be a power of two >= 1. It may be > 1 only when using a plane
// (flat surface) projection.
static const int ZOOM = 1;

static void CalculateConstants() {
    // Calculate the mesh dimensions in terms of pair of triangles.
    MESH_WIDTH = planet_elevation_map_texture->width / planet_elevation_map_detail_scale_down - 1;
    MESH_HEIGHT = planet_elevation_map_texture->height / planet_elevation_map_detail_scale_down - 1;
    // Calculate the number of elevation map samples used to calculate a vertex elevation
    // (which is guaranteed to leave no remainder).
    MESH_SAMPLE_AREA_WIDTH = planet_elevation_map_texture->width / (MESH_WIDTH + 1);
    MESH_SAMPLE_AREA_HEIGHT = planet_elevation_map_texture->height / (MESH_HEIGHT + 1);
    // Calculate the dimensions of the array of mesh models.
    WIDTH_IN_MESH_MODELS = (MESH_WIDTH + SUB_MESH_WIDTH - 1) / (SUB_MESH_WIDTH - 1);
    HEIGHT_IN_MESH_MODELS = (MESH_HEIGHT + SUB_MESH_HEIGHT - 1) / (SUB_MESH_HEIGHT - 1);
    // Calculate the offsets in terms of mesh vertices so that the mesh is centered
    // on (latitude, longitude). For a spherical planet model, these should be zero.
    if (planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
        X_OFFSET = 0;
        Y_OFFSET = 0;
    }
    else {
        X_OFFSET = (LONGITUDE + 180.0f) * (MESH_WIDTH + 1) / 360.0f - (MESH_WIDTH + 1) / 2 / ZOOM;
        Y_OFFSET = (LATITUDE + 90.0f) * (MESH_HEIGHT + 1) / 180.0f - (MESH_HEIGHT + 1) / 2 / ZOOM;
    }
}

static float CalculateTerrainDistance(float normalized_elevation_map_value) {
     // Calculate the distance from the planet center of the terrain vertex.
     return ((normalized_elevation_map_value -
         planet_elevation_scale->normalized_elevation_map_sea_level)
         * planet_elevation_scale->scale_factor * TERRAIN_HEIGHT_EXAGGERATION +
         planet_elevation_scale->sea_level_radius) * EXTRA_PLANET_SCALE_FACTOR;
}

// Create planet terrain using (WIDTH_IN_MESH_MODELS * SUBMESHES_Y) mesh models, using
// sreTexture::planet_elevation_map_texture with detail scale down factor
// elevation_map_detail_scale_down.

static void CalculatePlanetMesh(sreScene *scene, sreModel **mesh_model,
bool delete_texture_data) {
    // Create a vertex mesh by down scaling the elevation map.
    printf("Calculating vertices.\n");
    int v = 0;
    Point3D *vertex = new Point3D[MESH_WIDTH * MESH_HEIGHT + MESH_HEIGHT * 2 + 2];
    Point2D *texcoords = new Point2D[MESH_WIDTH * MESH_HEIGHT + MESH_HEIGHT * 2 + 2];
    int elevation_bits;
    float elevation_bits_normalization_factor;
    float planetary_radius = planet_elevation_scale->sea_level_radius;
    if (planet_elevation_map_texture->bit_depth == 16) {
        elevation_bits = 16;
        elevation_bits_normalization_factor = 1.0f / 65535.0f;
        // The 16-bit elevation map sea-level is at 186 for Earth.
//        elevation_offset = - 186.0f * Z_SCALE / 65535.0f;
    } 
    else {
        elevation_bits = 8;
        elevation_bits_normalization_factor = 1.0f / 255.0f;
        // The sea-level value for the 8-bit elevation maps is 23 for Earth.
//        elevation_offset = - 23.0f * Z_SCALE / 255.0f;
    }
    float X_SCALE = planetary_radius;
    float Y_SCALE = X_SCALE;
    for (int y = 0; y < MESH_HEIGHT; y++)
        for (int x = 0; x < MESH_WIDTH; x++) {
            // Calculate the average height in the area of size MESH_SAMPLE_AREA_WIDTH * MESH_SAMPLE_AREA_HEIGHT.
            double h = 0;
            for (int i = 0; i < MESH_SAMPLE_AREA_HEIGHT / ZOOM; i++)
                for (int j = 0; j < MESH_SAMPLE_AREA_WIDTH / ZOOM; j++) {
                    // Assume the first color component (red) is a value from 0 to 255 (or 65535)
                    // representing the height.
                    unsigned int elevation_value = planet_elevation_map_texture->LookupPixel(
                        x * (MESH_SAMPLE_AREA_WIDTH / ZOOM) + j + X_OFFSET * MESH_SAMPLE_AREA_WIDTH,
                        planet_elevation_map_texture->height - Y_OFFSET * MESH_SAMPLE_AREA_HEIGHT - (MESH_SAMPLE_AREA_HEIGHT / ZOOM) -
                        y * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + i)
                        & (((unsigned int)1 << elevation_bits) - 1);
                    h += elevation_value;
                }
            float longitude, latitude;
            float xcoord, ycoord, zcoord;
            // Calculate the elevation map value normalized to [0, 1].
            float normalized_elevation_value = elevation_bits_normalization_factor * h
                / ((MESH_SAMPLE_AREA_WIDTH / ZOOM) * (MESH_SAMPLE_AREA_HEIGHT / ZOOM));
            if (planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
                longitude = ((float)x + 0.5) / MESH_WIDTH * 2.0 * M_PI - M_PI;
                latitude = ((float)y + 0.5) / MESH_HEIGHT * M_PI - 0.5 * M_PI;
                float radius = CalculateTerrainDistance(normalized_elevation_value);
                xcoord = radius * cosf(latitude) * cosf(longitude);
                ycoord = radius * cosf(latitude) * sinf(longitude);
                zcoord = radius * sinf(latitude);
            }
            else {
               xcoord = (x - MESH_WIDTH / 2) * X_SCALE / MESH_WIDTH;
               ycoord = (y - MESH_HEIGHT / 2) * Y_SCALE / MESH_WIDTH;
               zcoord = (normalized_elevation_value -
                   planet_elevation_scale->normalized_elevation_map_sea_level)
                   * planet_elevation_scale->scale_factor * TERRAIN_HEIGHT_EXAGGERATION;
            }

            vertex[v].Set(xcoord, ycoord, zcoord);
            // Set the texcoords to the middle of the sampled area.
            // First calculate the texcoords in pixel units.
            float texcoords_x = x * (MESH_SAMPLE_AREA_WIDTH / ZOOM) + 0.5f * (MESH_SAMPLE_AREA_WIDTH / ZOOM) - 0.5 +
                X_OFFSET * MESH_SAMPLE_AREA_WIDTH;
            float texcoords_y = planet_elevation_map_texture->height - Y_OFFSET * MESH_SAMPLE_AREA_HEIGHT - (MESH_SAMPLE_AREA_HEIGHT / ZOOM) -
                y * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5f * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5;
            // The color/specularity/nightlight textures may have a different size, but the
            // normalized texture coordinates will be the same.
            texcoords[v].Set(
                texcoords_x / planet_elevation_map_texture->width,
                texcoords_y / planet_elevation_map_texture->height);
            v++;
        }

    int vertex_index_longitude_minus_180;
    int vertex_index_longitude_180;
    int vertex_index_latitude_minus_90;
    int vertex_index_latitude_90;
    if (!(planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL))
        goto skip_seams;
    // Define special column of vertices at - 180 degrees longitude.
    vertex_index_longitude_minus_180 = v;
    for (int y = 0; y < MESH_HEIGHT; y++) {
        vertex[v].Set(0.5 * (vertex[y * MESH_WIDTH].x + vertex[y * MESH_WIDTH + MESH_WIDTH - 1].x),
            0, 0.5 * (vertex[y * MESH_WIDTH].z + vertex[y * MESH_WIDTH + MESH_WIDTH - 1].z));
        float texcoords_y = planet_elevation_map_texture->height - Y_OFFSET * MESH_SAMPLE_AREA_HEIGHT - (MESH_SAMPLE_AREA_HEIGHT / ZOOM) -
            y * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5 * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5;
        texcoords[v].Set(0.0f,
            texcoords_y / planet_elevation_map_texture->height);
        v++;
    }
    // Define special column of vertices at 180 degrees longitude.
    vertex_index_longitude_180 = v;
    for (int y = 0; y < MESH_HEIGHT; y++) {
        vertex[v].Set(0.5 * (vertex[y * MESH_WIDTH].x + vertex[y * MESH_WIDTH + MESH_WIDTH - 1].x),
            0, 0.5 * (vertex[y * MESH_WIDTH].z + vertex[y * MESH_WIDTH + MESH_WIDTH - 1].z));
        float texcoords_y = planet_elevation_map_texture->height - Y_OFFSET * MESH_SAMPLE_AREA_HEIGHT - (MESH_SAMPLE_AREA_HEIGHT / ZOOM) -
            y * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5 * (MESH_SAMPLE_AREA_HEIGHT / ZOOM) + 0.5;
        texcoords[v].Set(1.0f,
            texcoords_y / planet_elevation_map_texture->height);
        v++;
    }
    // Define special vertices at the south and north polar caps.
    vertex_index_latitude_minus_90 = v;
    {
    float normalized_elevation_value = elevation_bits_normalization_factor *
        planet_elevation_map_texture->LookupPixel(
            planet_elevation_map_texture->width / 2, planet_elevation_map_texture->height - 1);
    float radius = CalculateTerrainDistance(normalized_elevation_value);
    vertex[v].Set(0, 0, - radius);
    texcoords[v].Set(0, 1.0);
    v++;
    vertex_index_latitude_90 = v;
    normalized_elevation_value = elevation_bits_normalization_factor *
        planet_elevation_map_texture->LookupPixel(planet_elevation_map_texture->width / 2, 0);
    radius = CalculateTerrainDistance(normalized_elevation_value);
    vertex[v].Set(0, 0, radius);
    texcoords[v].Set(0, 0);
    v++;
    }
skip_seams :

    // Free the elevation map texture data array when it is no longer required.
    if (delete_texture_data)
        planet_elevation_map_texture->ClearData();

    // Create mesh models.
    sreModelTriangle *triangle = new sreModelTriangle[(MESH_WIDTH - 1) * (MESH_HEIGHT - 1) * 2];
    int t = 0;
    for (int y = 0; y < MESH_HEIGHT - 1; y++)
        for (int x = 0; x < MESH_WIDTH - 1; x++) {
            triangle[t].AssignVertices(y * MESH_WIDTH + x, y * MESH_WIDTH + x + 1, (y + 1) *
                MESH_WIDTH + x + 1);
            triangle[t + 1].AssignVertices(y * MESH_WIDTH + x, (y + 1) * MESH_WIDTH + x + 1,
                (y + 1) * MESH_WIDTH + x);
            t += 2;
        }
    printf("Calculating normals.\n");
    for (int i = 0; i < t; i++) {
        triangle[i].normal = CalculateNormal(vertex[triangle[i].vertex_index[0]],
            vertex[triangle[i].vertex_index[1]], vertex[triangle[i].vertex_index[2]]);
        if (planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
            float sqmag = SquaredMag(triangle[i].normal);
            if (sqmag < 0.999 || sqmag > 1.001) {
                PlanetApplicationError("Invalid normal encountered.");
            }
        }
    }
    Vector3D *vertex_normal = new Vector3D[MESH_WIDTH * MESH_HEIGHT + MESH_HEIGHT * 2 + 2];
    for  (int y = 0; y < MESH_HEIGHT; y++)
        for (int x = 0; x < MESH_WIDTH; x++) {
            int t1 = y * (MESH_WIDTH - 1) * 2 + x * 2;
            int t2 = t1 + 1;
            Vector3D sum;
            sum.Set(0, 0, 0);
            if (x > 0) {
                if (y < MESH_HEIGHT - 1) {
                    sum += triangle[t1 - 2].normal;
                }
                if (y > 0) {
                   sum += triangle[t1 - (MESH_WIDTH - 1) * 2 - 2].normal;
                   sum += triangle[t2 - (MESH_WIDTH - 1) * 2 - 2].normal;
                }
            }
            if (x < MESH_WIDTH - 1) {
                if (y < MESH_HEIGHT - 1) {
                    sum += triangle[t1].normal;
                    sum += triangle[t2].normal;
                }
                if (y > 0)
                   sum += triangle[t2 - (MESH_WIDTH - 1) * 2].normal;
            }
            vertex_normal[y * MESH_WIDTH + x] = sum.Normalize();
            if (planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
                float sqmag = SquaredMag(vertex_normal[y * MESH_WIDTH + x]) ;
                if (sqmag < 0.999f || sqmag > 1.001f)
                   PlanetApplicationError("Invalid vertex normal encountered.");
            }
        }
    delete [] triangle;
    if (planet_mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
         goto skip_seam_normals;
    // Calculate vertex normals for special columns of vertices at - 180 and 180 degrees longitude.
    for (int y = 0; y < MESH_HEIGHT; y++) {
        vertex_normal[vertex_index_longitude_minus_180 + y] = (vertex_normal[y * MESH_WIDTH] +
            vertex_normal[y * MESH_WIDTH + MESH_WIDTH - 1]).Normalize();
        vertex_normal[vertex_index_longitude_180 + y] = (vertex_normal[y * MESH_WIDTH] +
            vertex_normal[y * MESH_WIDTH + MESH_WIDTH - 1]).Normalize();
    }
    vertex_normal[vertex_index_latitude_minus_90].Set(0, 0, - 1.0);
    vertex_normal[vertex_index_latitude_90].Set(0, 0, 1.0);
skip_seam_normals :

    printf("Assigning submeshes.\n");
    int total_triangle_count = 0;
    int total_triangle_count_reduced = 0;
    for (int sub_mesh_y = 0; sub_mesh_y < HEIGHT_IN_MESH_MODELS; sub_mesh_y++)
        for (int sub_mesh_x = 0; sub_mesh_x < WIDTH_IN_MESH_MODELS; sub_mesh_x++) {
            sreModel *model = mesh_model[MESH_INDEX(sub_mesh_x, sub_mesh_y)];
            sreLODModel *m = sreNewLODModel();
            int w = SUB_MESH_WIDTH;
            int h = SUB_MESH_HEIGHT;
            
            int x_offset = 0;
            int y_offset = 0;
#ifdef SPHERE
            // At longitude - 180 degrees, we need special vertices to cover the gap to 180 degrees.
            if (sub_mesh_x == 0) {
                w++;
                x_offset = 1;
            }
#endif
            if (sub_mesh_x * (SUB_MESH_WIDTH - 1) + w > MESH_WIDTH) {
                w = MESH_WIDTH - sub_mesh_x * (SUB_MESH_WIDTH - 1);
#ifdef SPHERE
                // Similarly, at longitude 180 degrees, we need to cover the gap to - 180 degrees.
                w++;
#endif
            }
#ifdef SPHERE
            if (sub_mesh_y == 0) {
                h++;
                y_offset = 1;
            }
#endif
            if (sub_mesh_y * (SUB_MESH_HEIGHT - 1) + h > MESH_HEIGHT) {
                h = MESH_HEIGHT - sub_mesh_y * (SUB_MESH_HEIGHT - 1);
#ifdef SPHERE
                // At latitude 90 degrees, we need to cover the gap to the north polar cap.
                h++;
#endif
            }
            // Use a smaller number of longitudinal vertices at higher latitudes.
            // This requires more work and is disabled.
            int longitudinal_step = 1;
#if defined(SPHERE) && false
            if (90.0f * sub_mesh_y * (SUB_MESH_HEIGHT - 1) / MESH_HEIGHT > 60.0f ||
            90.f * (sub_mesh_y * (SUB_MESH_HEIGHT - 1) + SUB_MESH_HEIGHT - 1) < - 60.0f)
                longitudual_step = 2;
#endif
            m->nu_vertices = w * h;
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
            m->nu_vertices += w * 2 + w * 2 + h * 2 + h * 2 + 8 + 4; // Extra edge vertices.
#endif
            m->vertex = dstNewAligned <Point3DPadded>(m->nu_vertices, 16);
            m->texcoords = new Point2D[m->nu_vertices];
            m->vertex_normal = new Vector3D[m->nu_vertices];
            float zmin = FLT_MAX;
            for (int y = 0; y < h; y++)
                for (int x = 0;; x += longitudinal_step) {
                    // Make sure the rightmost vertex is always included.
                    if (longitudinal_step > 1 && x >= w && x < w - 1 + longitudinal_step)
                        x = w - 1;
                    else if (x >= w)
                        break;
                    int mesh_x = sub_mesh_x * (SUB_MESH_WIDTH - 1);
                    int mesh_y = sub_mesh_y * (SUB_MESH_HEIGHT - 1);
                    int index;
#ifdef SPHERE
                    if (mesh_y + y == 0)
                        // South polar cap vertex.
                        index = vertex_index_latitude_minus_90;
                    else
                    if (mesh_y + y == MESH_HEIGHT)
                        // North polar cap vertex.
                        index = vertex_index_latitude_90;
                    // Special check for sphere to link up both sides at 180 / -180 degrees longitude.
                    else
                    if (mesh_x + x == 0)
                        // Use one of the extra vertices defined at longitude - 180 degrees.
                        index = vertex_index_longitude_minus_180 + (mesh_y + y);
                    else
                    if (mesh_x + x == MESH_WIDTH)
                        // Use one of the extra vertices defined at longitude 180 degrees.
                        index = vertex_index_longitude_180 + (mesh_y + y);
                    else
                        index = (mesh_y + y - y_offset) * MESH_WIDTH + (mesh_x + x - x_offset);
#else
                    index = (mesh_y + y) * MESH_WIDTH + (mesh_x + x);
#endif
                    if (index < 0 || index >= MESH_WIDTH * MESH_HEIGHT + MESH_HEIGHT * 2 + 2) {
                        printf("Index out of bounds (%d).\n", index);
                        exit(1);
                    }
                    m->vertex[y * w + x] = vertex[index];
                    m->texcoords[y * w + x] = texcoords[index];
                    m->vertex_normal[y * w + x] = vertex_normal[index];
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
                    if (vertex[index].z < zmin)
                        zmin = vertex[index].z;
#endif
                }
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
            zmin -= 0.001f * planetary_radius;
            if (zmin <= 0) {
                printf("zmin out of bounds.\n");
                exit(1);
            }
            // Duplicate vertices at the sides and add vertices on the ground.
            // Side at y = 0.
            int v = w * h;
            int *vertex_index_ymin_ground = new int[w];
            int *vertex_index_ymin_edge = new int[w];
            for (int x = 0; x < w; x++) {
                m->vertex[v] = m->vertex[x];
                m->vertex_normal[v] = Vector3D(0, -1.0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_ymin_edge[x] = v;
                v++;
                m->vertex[v] = Point3D(m->vertex[x].x, m->vertex[x].y, zmin);
                m->vertex_normal[v] = Vector3D(0, -1.0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_ymin_ground[x] = v;
                v++;
            }
            // Side at y = h - 1.
            int *vertex_index_ymax_ground = new int[w];
            int *vertex_index_ymax_edge = new int[w];
            for (int x = 0; x < w; x++) {
                m->vertex[v] = m->vertex[(h - 1) * w + x];
                m->vertex_normal[v] = Vector3D(0, 1.0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_ymax_edge[x] = v;
                v++;
                m->vertex[v] = Point3D(m->vertex[(h - 1) * w + x].x, m->vertex[(h - 1) * w + x].y, zmin);
                m->vertex_normal[v] = Vector3D(0, 1.0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_ymax_ground[x] = v;
                v++;
            }
            // Side at x = 0.
            int *vertex_index_xmin_ground = new int[h];
            int *vertex_index_xmin_edge = new int[h];
            for (int y = 0; y < h; y++) {
                m->vertex[v] = m->vertex[y * w];
                m->vertex_normal[v] = Vector3D(- 1.0, 0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_xmin_edge[y] = v;
                v++;
                m->vertex[v] = Point3D(m->vertex[y * w].x, m->vertex[y * w].y, zmin);
                m->vertex_normal[v] = Vector3D(- 1.0, 0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_xmin_ground[y] = v;
                v++;
            }
            // Side at x = w - 1.
            int *vertex_index_xmax_ground = new int[h];
            int *vertex_index_xmax_edge = new int[h];
            for (int y = 0; y < h; y++) {
                m->vertex[v] = m->vertex[y * w + w - 1];
                m->vertex_normal[v] = Vector3D(1.0, 0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_xmax_edge[y] = v;
                v++;
                m->vertex[v] = Point3D(m->vertex[y * w + w - 1].x, m->vertex[y * w + w - 1].y, zmin);
                m->vertex_normal[v] = Vector3D(1.0, 0, 0);
                m->texcoords[v] = Point2D(0, 0);
                vertex_index_xmax_ground[y] = v;
                v++;
            }
            // Create bottom corner vertices for each side.
            int ymin_corner_vertex_index = v;
            m->vertex[v] =  Point3D(m->vertex[vertex_index_ymin_ground[0]].x, m->vertex[vertex_index_ymin_ground[0]].y, 0);
            m->vertex_normal[v] = Vector3D(0, - 1.0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] = Point3D(m->vertex[vertex_index_ymin_ground[w - 1]].x,
                m->vertex[vertex_index_ymin_ground[w - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(0, - 1.0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            int ymax_corner_vertex_index = v;
            m->vertex[v] = Point3D(m->vertex[vertex_index_ymax_ground[w - 1]].x,
                m->vertex[vertex_index_xmax_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 1.0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] =  Point3D(m->vertex[vertex_index_ymax_ground[0]].x,
                m->vertex[vertex_index_xmin_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 1.0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            int xmin_corner_vertex_index = v;             
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmin_ground[h - 1]].x,
                m->vertex[vertex_index_xmin_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(- 1.0, 0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] =  Point3D(m->vertex[vertex_index_xmin_ground[0]].x, m->vertex[vertex_index_xmin_ground[0]].y, 0);
            m->vertex_normal[v] = Vector3D(- 1.0, 0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            int xmax_corner_vertex_index = v;
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmax_ground[0]].x, m->vertex[vertex_index_xmax_ground[0]].y, 0);
            m->vertex_normal[v] = Vector3D(1.0, 0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmax_ground[h - 1]].x,
                m->vertex[vertex_index_xmax_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(1.0, 0, 0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            // Create bottom corner vertices for bottom.
            int bottom_corner_vertex_index = v;
            m->vertex[v] =  Point3D(m->vertex[vertex_index_xmin_ground[0]].x, m->vertex[vertex_index_xmin_ground[0]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 0, - 1.0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmax_ground[0]].x, m->vertex[vertex_index_xmax_ground[0]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 0, - 1.0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmin_ground[h - 1]].x,
                m->vertex[vertex_index_xmin_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 0, - 1.0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            m->vertex[v] = Point3D(m->vertex[vertex_index_xmax_ground[h - 1]].x,
                m->vertex[vertex_index_xmax_ground[h - 1]].y, 0);
            m->vertex_normal[v] = Vector3D(0, 0, - 1.0);
            m->texcoords[v] = Vector2D(0, 0);
            v++;
            if (v != m->nu_vertices) {
                printf("Mismatch in the number of vertices (%d vs %d).\n", v, m->nu_vertices);
                exit(1);
            }
#endif
            m->nu_triangles = 2 * (w - 1) * (h - 1);
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
            // Extra triangles at the y = 0 side.
            int side_triangles_ymin = 1 + (w - 3) / 2 + ((w - 3) / 2) * 2 + 3;
            if ((w - 3) & 1)
                side_triangles_ymin++; // One extra triangle in the middle section (x == uneven).
            m->nu_triangles += side_triangles_ymin;
            // Triangle fan. Number of vertices in the triangle fan is equal to the number of side vertices created + 1.
            // Side vertices are created at the start (two), one when x is even.
            // is uneven.
            m->nu_triangles += 2 + (w - 3) / 2 + 1;
            // Extra triangles at the y = h - 1 side.
            int side_triangles_ymax = 2 + (w - 3) / 2 + ((w - 3) / 2) * 2 + 1;
            if ((w - 3) & 1)
                side_triangles_ymax += 2; // Two extra triangles in the middle section (x == even)
            if (!((w - 3) & 1))
                side_triangles_ymax++; // One extra triangle at the end.
            m->nu_triangles += side_triangles_ymax;
            // Triangle fan. Number of vertices in the triangle fan is equal to the number of side vertices created + 1.
            // Side vertices are created at the start (two), one when x is even.
            m->nu_triangles += 2 + (w - 3) / 2 + 1;
            // Extra triangles at the x = 0 side.
            m->nu_triangles += 2 + (h - 3) / 2 + ((h - 3) / 2) * 2 + 1;
            if ((h - 3) & 1)
                m->nu_triangles += 2; // Two extra triangles in the middle section (y == even)
            if (!((h - 3) & 1))
                m->nu_triangles++; // One extra triangle at the end.
            // Triangle fan.
            m->nu_triangles += 2 + (h - 3) / 2 + 1;
            // Extra triangles at the x = w - 1 side.
            m->nu_triangles += 1 + (h - 3) / 2 + ((h - 3) / 2) * 2 + 3;
            if ((h - 3) & 1)
                m->nu_triangles++; // One extra triangle in the middle section (y == uneven)
            // Triangle fan.
            m->nu_triangles += 2 + (h - 3) / 2 + 1;
            m->nu_triangles += 2; // Two triangles for the bottom.
            bool *ymin_ground_vertex_used = new bool[w];
            bool *ymax_ground_vertex_used = new bool[w];
            bool *xmin_ground_vertex_used = new bool[h];
            bool *xmax_ground_vertex_used = new bool[h];
            for (int i = 0; i < w; i++) {
                ymin_ground_vertex_used[i] = false;
                ymax_ground_vertex_used[i] = false;
            }
            for (int i = 0; i < h; i++) {
                xmin_ground_vertex_used[i] = false;
                xmax_ground_vertex_used[i] = false;
            }
            // Create side triangles.
            int counted_triangles_ymin = 0;
            int counted_triangles_ymax = 0;
#endif
            m->triangle = new sreModelTriangle[m->nu_triangles];
            int t = 0;
            for (int y = 0; y < h - 1; y++)
                for (int x = 0; x < w - 1; x++) {
                    m->triangle[t].AssignVertices(y * w + x, y * w + x + 1, (y + 1) * w + x + 1);
                    m->triangle[t + 1].AssignVertices(y * w + x, (y + 1) * w + x + 1, (y + 1) * w + x);
                    t += 2;
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
                    // If we are on a side, add triangles to make a closed volume for shadow volumes.
                    if (y == 0) {
                        if (x == 0) {
                            m->triangle[t].AssignVertices(vertex_index_ymin_edge[x], vertex_index_ymin_ground[x],
                                vertex_index_ymin_ground[x + 1]);
                            ymin_ground_vertex_used[x] = true;
                            ymin_ground_vertex_used[x + 1] = true;
                            t++;
                            counted_triangles_ymin++;
                        }
                        else
                        if (x & 1) {
                            m->triangle[t].AssignVertices(vertex_index_ymin_edge[x], vertex_index_ymin_edge[x - 1],
                                vertex_index_ymin_ground[x]);
                            ymin_ground_vertex_used[x] = true;
                            t++;
                            counted_triangles_ymin++;
                            if (x == w - 2) {
                                m->triangle[t].AssignVertices(vertex_index_ymin_edge[x + 1], vertex_index_ymin_edge[x],
                                    vertex_index_ymin_ground[x]);
                                m->triangle[t + 1].AssignVertices(vertex_index_ymin_edge[x + 1], vertex_index_ymin_ground[x],
                                    vertex_index_ymin_ground[x + 1]);
                                ymin_ground_vertex_used[x + 1] = true;
                                t += 2;
                                counted_triangles_ymin += 2;
                            }
                        }
                        else {
                            m->triangle[t].AssignVertices(vertex_index_ymin_edge[x], vertex_index_ymin_edge[x - 1],
                                vertex_index_ymin_ground[x - 1]);
                            m->triangle[t + 1].AssignVertices(vertex_index_ymin_edge[x], vertex_index_ymin_ground[x - 1],
                                vertex_index_ymin_ground[x + 1]);
                            ymin_ground_vertex_used[x - 1] = true;
                            ymin_ground_vertex_used[x + 1] = true;
                            t += 2;
                            counted_triangles_ymin += 2;
                            if (x == w - 2) {
                                m->triangle[t].AssignVertices(vertex_index_ymin_edge[x + 1], vertex_index_ymin_edge[x],
                                    vertex_index_ymin_ground[x + 1]);
                                t++;
                                counted_triangles_ymin++;
                            }
                        }
                    }
                    if (y == h - 2) {
                        if (x == w - 2) {
                            m->triangle[t].AssignVertices(vertex_index_ymax_edge[x + 1], vertex_index_ymax_ground[x + 1],
                                vertex_index_ymax_ground[x]);
                            m->triangle[t + 1].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_edge[x + 1],
                                vertex_index_ymax_ground[x]);
                            ymax_ground_vertex_used[x] = true;
                            ymax_ground_vertex_used[x + 1] = true;
                            t += 2;
                            counted_triangles_ymax += 2;
                        }
                        else
                        if (x == 0 && ((w - 3) & 1)) {
                            m->triangle[t].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_edge[x + 1],
                                vertex_index_ymax_ground[x]);
                            ymax_ground_vertex_used[x] = true;
                            t++;
                            counted_triangles_ymax++;
                        }
                        else
                        if (x == 0 && !((w - 3) & 1)) {
                            m->triangle[t].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_edge[x + 1],
                                vertex_index_ymax_ground[x]);
                            m->triangle[t + 1].AssignVertices(vertex_index_ymax_edge[x + 1], vertex_index_ymax_ground[x + 1],
                                vertex_index_ymax_ground[x]);
                            ymax_ground_vertex_used[x] = true;
                            ymax_ground_vertex_used[x + 1] = true;
                            t += 2;
                            counted_triangles_ymax += 2;
                        }
                        else
                        if ((w - 3 - x) & 1) {
                            m->triangle[t].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_edge[x + 1],
                                vertex_index_ymax_ground[x]);
                            ymax_ground_vertex_used[x] = true;
                            t++;
                            counted_triangles_ymax++;
                        }
                        else {
                            m->triangle[t].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_edge[x + 1],
                                vertex_index_ymax_ground[x + 1]);
                            m->triangle[t + 1].AssignVertices(vertex_index_ymax_edge[x], vertex_index_ymax_ground[x + 1],
                                vertex_index_ymax_ground[x - 1]);
                            ymax_ground_vertex_used[x - 1] = true;
                            ymax_ground_vertex_used[x + 1] = true;
                            t += 2;
                            counted_triangles_ymax += 2;
                        }
                    }
                    if (x == 0) {
                        if (y == h - 2) {
                            m->triangle[t].AssignVertices(vertex_index_xmin_edge[y + 1], vertex_index_xmin_ground[y + 1],
                                vertex_index_xmin_ground[y]);
                            m->triangle[t + 1].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_edge[y + 1],
                                vertex_index_xmin_ground[y]);
                            xmin_ground_vertex_used[y] = true;
                            xmin_ground_vertex_used[y + 1] = true;
                            t += 2;
                        }
                        else
                        if (y == 0 && ((h - 3) & 1)) {
                            m->triangle[t].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_edge[y + 1],
                                vertex_index_xmin_ground[y]);
                            xmin_ground_vertex_used[y] = true;
                            t++;
                        }
                        else
                        if (y == 0 && !((h - 3) & 1)) {
                            m->triangle[t].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_edge[y + 1],
                                vertex_index_xmin_ground[y]);
                            m->triangle[t + 1].AssignVertices(vertex_index_xmin_edge[y + 1], vertex_index_xmin_ground[y + 1],
                                vertex_index_xmin_ground[y]);
                            xmin_ground_vertex_used[y] = true;
                            xmin_ground_vertex_used[y + 1] = true;
                            t += 2;
                        }
                        else
                        if ((h - 3 - y) & 1) {
                            m->triangle[t].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_edge[y + 1],
                                vertex_index_xmin_ground[y]);
                            xmin_ground_vertex_used[y] = true;
                            t++;
                        }
                        else {
                            m->triangle[t].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_edge[y + 1],
                                vertex_index_xmin_ground[y + 1]);
                            m->triangle[t + 1].AssignVertices(vertex_index_xmin_edge[y], vertex_index_xmin_ground[y + 1],
                                vertex_index_xmin_ground[y - 1]);
                            xmin_ground_vertex_used[y - 1] = true;
                            xmin_ground_vertex_used[y + 1] = true;
                            t += 2;
                        }
                    }
                    if (x == w - 2) {
                        if (y == 0) {
                            m->triangle[t].AssignVertices(vertex_index_xmax_edge[y], vertex_index_xmax_ground[y],
                                vertex_index_xmax_ground[y + 1]);
                            xmax_ground_vertex_used[y] = true;
                            xmax_ground_vertex_used[y + 1] = true;
                            t++;
                        }
                        else
                        if (y & 1) {
                            m->triangle[t].AssignVertices(vertex_index_xmax_edge[y], vertex_index_xmax_edge[y - 1],
                                vertex_index_xmax_ground[y]);
                            xmax_ground_vertex_used[y] = true;
                            t++;
                            if (y == h - 2) {
                                m->triangle[t].AssignVertices(vertex_index_xmax_edge[y + 1], vertex_index_xmax_edge[y],
                                    vertex_index_xmax_ground[y]);
                                m->triangle[t + 1].AssignVertices(vertex_index_xmax_edge[y + 1], vertex_index_xmax_ground[y],
                                    vertex_index_xmax_ground[y + 1]);
                                xmax_ground_vertex_used[y + 1] = true;
                                t += 2;
                            }
                        }
                        else {
                            m->triangle[t].AssignVertices(vertex_index_xmax_edge[y], vertex_index_xmax_edge[y - 1],
                                vertex_index_xmax_ground[y - 1]);
                            m->triangle[t + 1].AssignVertices(vertex_index_xmax_edge[y], vertex_index_xmax_ground[y - 1],
                                vertex_index_xmax_ground[y + 1]);
                            xmax_ground_vertex_used[y - 1] = true;
                            xmax_ground_vertex_used[y + 1] = true;
                            t += 2;
                            if (y == h - 2) {
                                m->triangle[t].AssignVertices(vertex_index_xmax_edge[y + 1], vertex_index_xmax_edge[y],
                                    vertex_index_xmax_ground[y + 1]);
                                t++;
                            }
                        }
                    }
#endif
                }
#ifdef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
            if (counted_triangles_ymin != side_triangles_ymin) {
                printf("Mismatch in number of ymin side triangles (%d vs %d).\n", counted_triangles_ymin,
                    side_triangles_ymin);
                exit(1);
            }
            if (counted_triangles_ymax != side_triangles_ymax) {
                printf("Mismatch in number of ymax side triangles (%d vs %d).\n", counted_triangles_ymax,
                    side_triangles_ymax);
                exit(1);
            }
            // Create triangle fans for the bottom part of each side.
            int previous_vertex = vertex_index_ymin_ground[0];
            for (int x = 1; x < w; x++)
                if (ymin_ground_vertex_used[x]) {
                    m->triangle[t].AssignVertices(ymin_corner_vertex_index, vertex_index_ymin_ground[x], previous_vertex);
                    previous_vertex = vertex_index_ymin_ground[x];
                    t++;
                }
            m->triangle[t].AssignVertices(ymin_corner_vertex_index, ymin_corner_vertex_index + 1, previous_vertex);
            t++;
            previous_vertex = vertex_index_ymax_ground[w - 1];
            for (int x = w - 2; x >= 0; x--)
                if (ymax_ground_vertex_used[x]) {
                    m->triangle[t].AssignVertices(ymax_corner_vertex_index, vertex_index_ymax_ground[x], previous_vertex);
                    previous_vertex = vertex_index_ymax_ground[x];
                    t++;
                }
            m->triangle[t].AssignVertices(ymax_corner_vertex_index, ymax_corner_vertex_index + 1, previous_vertex);
            t++;
            previous_vertex = vertex_index_xmin_ground[h - 1];
            for (int y = h - 2; y >= 0; y--)
                if (xmin_ground_vertex_used[y]) {
                    m->triangle[t].AssignVertices(xmin_corner_vertex_index, vertex_index_xmin_ground[y], previous_vertex);
                    previous_vertex = vertex_index_xmin_ground[y];
                    t++;
                }
            m->triangle[t].AssignVertices(xmin_corner_vertex_index, xmin_corner_vertex_index + 1, previous_vertex);
            t++;
            previous_vertex = vertex_index_xmax_ground[0];
            for (int y = 1; y < h; y++)
                if (xmax_ground_vertex_used[y]) {
                    m->triangle[t].AssignVertices(xmax_corner_vertex_index, vertex_index_xmax_ground[y], previous_vertex);
                    previous_vertex = vertex_index_xmax_ground[y];
                    t++;
                }
            m->triangle[t].AssignVertices(xmax_corner_vertex_index, xmax_corner_vertex_index + 1, previous_vertex);
            t++;
            // Create the bottom.
            m->triangle[t].AssignVertices(bottom_corner_vertex_index, bottom_corner_vertex_index + 2,
                bottom_corner_vertex_index + 1);
            m->triangle[t + 1].AssignVertices(bottom_corner_vertex_index + 2, bottom_corner_vertex_index + 3,
                bottom_corner_vertex_index + 1);
            t += 2;
            if (t != m->nu_triangles) {
                printf("Mismatch in the number of triangles (%d vs %d), w = %d, h = %d.\n", t, m->nu_triangles, w, h);
                exit(1);
            }
            delete [] vertex_index_ymin_ground;
            delete [] vertex_index_ymin_edge;
            delete [] vertex_index_ymax_ground;
            delete [] vertex_index_ymax_edge;
            delete [] vertex_index_xmin_ground;
            delete [] vertex_index_xmin_edge;
            delete [] vertex_index_xmax_ground;
            delete [] vertex_index_xmax_edge;
            delete [] ymin_ground_vertex_used;
            delete [] ymax_ground_vertex_used;
            delete [] xmin_ground_vertex_used;
            delete [] xmax_ground_vertex_used;
#endif
            // When the vertex colors flag is set, do a lookup of colors from the texture and
            // delete texcoords.
            if (planet_mesh_flags & MESH_FLAG_VERTEX_COLORS) {
                m->colors = new Color[m->nu_vertices];
                for (int i = 0; i < m->nu_vertices; i++) {
                    float u = m->texcoords[i].x;
                    float v = m->texcoords[i].y;
//                    sreMessage(SRE_MESSAGE_INFO, "u = %f, v = %f", u, v);
                    Color c;
                    planet_color_map_texture->TextureLookupNearest(u, v, c);
                    m->colors[i] = c;
                }
                delete [] m->texcoords;
            }
            else {
                // Otherwise, flip the texcoords y coordinates (required so that latitudes are
                // not flipped).
                for (int i = 0; i < m->nu_vertices; i++)
                    m->texcoords[i].y = 1.0f - m->texcoords[i].y;
            }

            m->flags = SRE_POSITION_MASK | SRE_NORMAL_MASK;
            if (planet_mesh_flags & MESH_FLAG_VERTEX_COLORS)
                m->flags |= SRE_COLOR_MASK;
            else
                m->flags |= SRE_TEXCOORDS_MASK;
#ifndef CLOSED_SEGMENTS_FOR_SHADOW_VOLUMES
//            m->flags |= SRE_LOD_MODEL_NO_SHADOW_VOLUME_SUPPORT;
            m->flags |= SRE_LOD_MODEL_NOT_CLOSED;
#endif
            m->RemoveEmptyTriangles();
            m->RemoveUnusedVertices();
            m->CalculateTriangleNormals();
            total_triangle_count += m->nu_triangles;
            m->ReduceTriangleCount(0.5f, 0.05f, true, 0.995); // Cost threshold was 0.03.
            total_triangle_count_reduced += m->nu_triangles;
            m->SortVerticesOptimalDimension();
            // Vertex normals cannot be recalculated because it would result in discrepancies at the edges.
            m->CalculateTriangleNormals();
            model->lod_model[0] = m;
            model->lod_model[1] = sreNewLODModel();
            m->Clone(model->lod_model[1]);
            model->lod_model[1]->ReduceTriangleCount(0.2f, 0.4f, true, 0.97f);
            model->lod_model[1]->SortVerticesOptimalDimension();
            model->lod_model[1]->CalculateTriangleNormals();
            // Create a third LOD model.
            model->lod_model[2] = sreNewLODModel();
            model->lod_model[1]->Clone(model->lod_model[2]);
            model->lod_model[2]->ReduceTriangleCount(0.1f, 0.8f, true, 0.9f);
            model->lod_model[2]->SortVerticesOptimalDimension();
            model->lod_model[2]->CalculateTriangleNormals();
            model->nu_lod_levels = 3;

            // Mandate a significant (> 30%) reduction in triangle count between
            // LOD levels.
            float ratio2_0 = (float)model->lod_model[2]->nu_triangles /
                (float)model->lod_model[0]->nu_triangles;
            float ratio1_0 = (float)model->lod_model[1]->nu_triangles /
                (float)model->lod_model[0]->nu_triangles;
            float ratio2_1 = (float)model->lod_model[2]->nu_triangles /
                (float)model->lod_model[1]->nu_triangles;
            model->lod_threshold_scaling = 1.0f;
            if (ratio2_0 >= 0.7f) {
                model->nu_lod_levels = 1;
                delete model->lod_model[1];
                delete model->lod_model[2];
            }
            else if (ratio1_0 < 0.7f) {
                if (ratio2_1 < 0.7f) {
                    // Both LOD levels reduce the count noticeably, keep all three.
                }
                else {
                    // Discard level 2 since it does not offer a significantly reduced
                    // count.
                    delete model->lod_model[2];
                    model->nu_lod_levels = 2;
                }
            }
            else {
                // Only the ratio between levels 2 and 0 reaches 70%.
                // Use level 2 as level 1, and increase threshold scaling so that
                // level 1 is triggered as would otherwise be level 2.
                delete model->lod_model[1];
                model->lod_model[1] = model->lod_model[2];
                model->nu_lod_levels = 2;
                model->lod_threshold_scaling =
                    SRE_LOD_LEVEL_1_THRESHOLD / SRE_LOD_LEVEL_2_THRESHOLD;
            }
            printf("Using %d out of 3 LOD levels.\n", model->nu_lod_levels);

            for (int i = 0; i < model->nu_lod_levels; i++)
                model->lod_model[i]->cache_coherency_sorting_hint = SRE_SORTING_HINT_DO_NOT_SORT;

            model->CalculateBounds();
            model->collision_shape_static = SRE_COLLISION_SHAPE_STATIC;
            model->collision_shape_dynamic = SRE_COLLISION_SHAPE_CONVEX_HULL;
            scene->RegisterModel(model);
//            printf("%d vertices, %d triangles.\n", m->nu_vertices, m->nu_triangles);
        }
    delete [] vertex;
    delete [] texcoords;
    delete [] vertex_normal;
    printf("%d of %d triangles (%d%%) removed by edge collapse.\n", total_triangle_count - total_triangle_count_reduced,
        total_triangle_count, (total_triangle_count - total_triangle_count_reduced) * 100 / total_triangle_count);
    return;
}

// Loading/saving of Planet model meshes.

static char *GetMeshModelFileName(int x, int y) {
    char *filename = new char[256];
    char *stripped_filename = strchr((char *)planet_elevation_map_filename, '/');
    if (stripped_filename != NULL)
        stripped_filename++;
    else
        stripped_filename = (char *)planet_elevation_map_filename;
    char *attribute_type_str;
    if (planet_mesh_flags & MESH_FLAG_VERTEX_COLORS)
        attribute_type_str = "colors";
    else
        attribute_type_str = "texcoords";
    sprintf(filename, "planet-meshes/planet-mesh-x%dy%d-elevation-map-%s-%dx%d-detail-%d-attribute-%s"
        ".srebinarymodel",
        x, y, stripped_filename,
        planet_elevation_map_texture->width, planet_elevation_map_texture->height,
        planet_elevation_map_detail_scale_down,
        attribute_type_str);
    char *path = new char[256];
    int path_size = 256;
    while (getcwd(path, path_size) == NULL) {
        delete [] path;
        path_size *= 2;
        path = new char[path_size];
    }
    char *full_path = new char[strlen(path) + strlen(filename) + 2];
    strcpy(full_path, path);
    delete [] path;
    strcat(full_path, "/");
    strcat(full_path, filename);
    delete [] filename;
    return full_path;
}

static bool FileExists(const char *filename) {
#ifdef __GNUC__
    struct stat stat_buf;
    int r = stat(filename, &stat_buf);
    if (r == - 1)
        return false;
    return true;
#else
    FILE *f = fopen(filename, "rb");
    if (f == NULL)
	return false;
    fclose(f);
    return true;
#endif
}

static bool MeshObjectFilesExist() {
    for (int y = 0; y < HEIGHT_IN_MESH_MODELS; y++)
        for (int x = 0; x < WIDTH_IN_MESH_MODELS; x++) {
            char *filename = GetMeshModelFileName(x, y);
            bool exists = FileExists(filename);
            delete [] filename;
            if (!exists)
                return false;
        }
    return true;
}

static bool LoadMeshObjects(sreScene *scene, sreModel **mesh_model) {
    if (!MeshObjectFilesExist())
        return false;
    for (int y = 0; y < HEIGHT_IN_MESH_MODELS; y++)
        for (int x = 0; x < WIDTH_IN_MESH_MODELS; x++) {
            char *filename = GetMeshModelFileName(x, y);
            mesh_model[MESH_INDEX(x, y)] =
                sreReadModelFromSREBinaryModelFile(scene, filename, 0);
            delete [] filename;
        }
    return true;
}

static void SaveMeshObjects(sreScene *scene, sreModel **mesh_model) {
    if (!FileExists("planet-meshes")) {
        system("mkdir planet-meshes");
    }
    for (int y = 0; y < HEIGHT_IN_MESH_MODELS; y++)
        for (int x = 0; x < WIDTH_IN_MESH_MODELS; x++) {
            char *filename = GetMeshModelFileName(x, y);
            sreSaveModelToSREBinaryModelFile(mesh_model[MESH_INDEX(x, y)], filename, 0);
            delete [] filename;
        }
}

// Create a planet mesh based on elevation_map_texture. The texture and its data storage
// are freed during processing to limit memory usage.

void CreatePlanetMesh(sreScene *scene, const char *elevation_map_filename,
sreTexture *elevation_map_texture, int detail_division_factor,
sreTexture *color_map_texture, int _color_map_horizontal_subdivisions,
int _color_map_vertical_subdivisions, unsigned int mesh_flags, int& nu_models,
sreModel **&mesh_model) {
    planet_elevation_map_filename = elevation_map_filename;
    planet_elevation_map_texture = elevation_map_texture;
    planet_elevation_map_detail_scale_down = detail_division_factor;
    planet_color_map_texture = color_map_texture;
    color_map_horizontal_subdivisions = _color_map_horizontal_subdivisions;
    color_map_vertical_subdivisions = _color_map_vertical_subdivisions;
    planet_mesh_flags = mesh_flags;

    CalculateConstants();
    nu_models = WIDTH_IN_MESH_MODELS * HEIGHT_IN_MESH_MODELS;
    mesh_model = new sreModel *[WIDTH_IN_MESH_MODELS * HEIGHT_IN_MESH_MODELS];
    planet_mesh_model = mesh_model;

    for (int y = 0; y < HEIGHT_IN_MESH_MODELS; y++)
        for (int x = 0; x < WIDTH_IN_MESH_MODELS; x++)
             mesh_model[MESH_INDEX(x, y)] = new sreModel;
    bool loaded = false;
    if (!(planet_mesh_flags & (MESH_FLAG_NO_CACHE | MESH_FLAG_NEW_CACHE)))
        loaded = LoadMeshObjects(scene, mesh_model);
    if (!loaded) {
        CalculatePlanetMesh(scene, mesh_model, true);
        if (!(planet_mesh_flags & MESH_FLAG_NO_CACHE))
            SaveMeshObjects(scene, mesh_model);
    }
}
