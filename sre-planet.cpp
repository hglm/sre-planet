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
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdarg.h>
#include <signal.h>

#include "sre.h"
#include "sreBackend.h"

#include "planet-config.h"
#include "sre-planet.h"

// MatrixDouble3D not yet implemented.
typedef Matrix3D MatrixDouble3D;

class PlanetApplication : public sreBulletPhysicsApplication {
public : 
    virtual void StepBeforeRender(double demo_time);
    virtual void StepBeforePhysics(double demo_time);
};

// Great circle class.

class GreatCircleRouteSpec {
public :
    double long0, lat0;
    double long1, lat1;
    float height0, height1;
    float pitch0, pitch1;
};

class GreatCircleRoute : public GreatCircleRouteSpec {
private :
    double long_diff;
    double azimuth[2];
    double central_angle;
    double nodal_azimuth;
    double nodal_angle[2];
    double nodal_longitude;

public :
    void CalculateGreatCircle();
    void CalculateIntermediatePoint(double t, double& longitude,
        double& latitude, double& azimuth);
};

void GreatCircleRoute::CalculateGreatCircle() {
    long_diff = (long1 - long0) * M_PI / 180.0;
    azimuth[0] = atan2(sin(long_diff), cos(lat0 * M_PI / 180.0) *
        tan(lat1 * M_PI / 180.0) - sin(lat0 * M_PI / 180.0) * cos(long_diff));
    central_angle = acos(sin(lat0 * M_PI / 180.0) * sin(lat1 * M_PI / 180.0) +
        cos(lat0 * M_PI / 180.0) * cos(lat1 * M_PI / 180.0)
        * cos(long_diff));
    nodal_azimuth = atan2(sin(azimuth[0]) * cos(lat0 * M_PI / 180.0),
        sqrt(pow(cos(azimuth[0]), 2.0) + pow(sin(azimuth[0]), 2.0) *
        pow(sin(lat0 * M_PI / 180.0), 2.0)));
    if (lat0 == 0 && azimuth[0] == 0.5 * M_PI)
        nodal_angle[0] = 0;
    else
        nodal_angle[0] = atan2(tan(lat0 * M_PI / 180.0), cos(azimuth[0]));
    nodal_angle[1] = nodal_angle[0] + central_angle;
    double lon_A_P0 = atan2(sin(nodal_azimuth) * sin(nodal_angle[0]), cos(nodal_angle[0]));
    nodal_longitude = long0 * M_PI / 180.0 - lon_A_P0;
}

void GreatCircleRoute::CalculateIntermediatePoint(double t, double& longitude,
double& latitude, double& azimuth) {
    double ang = nodal_angle[0] + t * (nodal_angle[1] - nodal_angle[0]); 
    latitude = atan2(cos(nodal_azimuth) * sin(ang),
        sqrt(pow(cos(ang), 2.0) + pow(sin(nodal_azimuth), 2.0) *
        pow(sin(ang), 2.0)));
    longitude = atan2(sin(nodal_azimuth) * sin(ang), cos(ang)) + nodal_longitude;
    azimuth = atan2(tan(nodal_azimuth), cos(ang));
}


static const char *elevation_map_filename;
static const char *color_map_filename;
static const char *specularity_map_filename;
static const char *night_map_filename;
static sreTexture *elevation_map_texture;
static sreTexture *color_map_texture;
static sreTexture *specularity_map_texture;
static sreTexture *night_map_texture;
static int texture_size;
static int elevation_map_detail_scale_down;
static unsigned int mesh_flags;
static bool night_only;
static bool no_text_overlay;
static bool seasonal_bmng;
static int planet_preset;

static sreScene *scene;
static int player_object = -1;
static int spacecraft_object = -1;
static int sun_object;
static int directional_light;
static int spacecraft_spot_light;
static float saved_hovering_height;
static MatrixDouble3D saved_spacecraft_rotation_matrix;

static float day_interval = DEFAULT_DAY_INTERVAL;
static bool display_time = true;
static bool physics = true;
static bool create_spacecraft = true;
#ifdef SHOW_SPACECRAFT
static bool show_spacecraft = true;
#else
static bool show_spacecraft = false;
#endif
static float sun_light_factor = 1.0f;
static float extra_lod_threshold_scaling = 5.0f;
static sreTexture *spacecraft_emission_map;

static int current_month;
static sreTexture *seasonal_color_map_texture;
static sreTexture *bmng_seasonal_color_map_texture[12];

// Planet data.

enum {
    PLANET_EARTH = 0,
    PLANET_MOON = 1,
    PLANET_MARS = 2,
    PLANET_MERCURY = 3,
    PLANET_EARTH_BMNG = 4
};

static const char *planet_preset_filename[][5] = {
    { "16_bit_elevation_10800x5400", "6_merged_color_ice_16384",
        "water_16384", "cities_16384" },
    { "moon-elevation-shifted-23040", "moon-color-monochrome-16384",
        NULL, NULL },
    { "mars_elevation_4096.png", NULL, NULL, NULL },
    { "mercury_elevation-23040x11520", "mercury-interpolated-reflectance",
        NULL, NULL },
    // Blue Marble Next Generation data set (month May). Visible Earth data set, NASA
    { "BMNG/srtm_ramp2.world.16384x8192", "BMNG/world.200405.3x21600x10800-16384",
        NULL, NULL }
};

// Elevation map scale factor, "Sea" level radius to planet center, and correction to
// elevation map value (assuming range [0, 1]) so that 0 maps to sea level.

#define REAL_EARTH_RADIUS 6378.1f
#define REAL_MOON_RADIUS 1737.4f
#define REAL_MARS_RADIUS 3396.2f
#define REAL_MERCURY_RADIUS 2440.0f

const ElevationScale planet_elevation_scale[] = {
    // Earth. Sea level at elevation map value 186. Elevation map is normalized to
    // [0, 65535.0f].
    { 10000.0f, (REAL_EARTH_RADIUS / 10000.0f) * 8.800f / (65535.0f - 186.0f) * 65535.0f,
      186 / 65535.0f },
    // Moon. Elevation map integer values represent 0.5m.
    { 10000.0f, (REAL_MOON_RADIUS / 10000.0f) * 0.0005f * 65535.0f,  18250.0f / 65535.0f },
    // Mars.
    { 10000.0f, (REAL_MARS_RADIUS / 10000.0f) * 0.0005f * 65535.0f,  18250.0f / 65535.0f },
    // Mercury. Elevation map in range [0, 39796], reference level 11778.
    { 10000.0f, (REAL_MERCURY_RADIUS / 10000.0f) * 10.4f / 39796 * 65535.0f,
        11778 / 65535.0f },
    // Earth, Blue Marble Next Generation data set (Visible Earth, NASA). Elevation map
    // is normalized to [0, 65535.0f].
    { 10000.0f, (REAL_EARTH_RADIUS / 10000.0f) * 8.800f / (65535.0f - 186.0f) * 65535.0f,
      186 / 65535.0f }
};

#define PLANETARY_RADIUS (planet_elevation_scale[planet_preset].sea_level_radius)

static void Fill(sreTexture *texture, int x, int y, int w, int h, unsigned int color) {
    for (int i = y; i < y + h; i++)
        for (int j = x; j < x + w; j++)
            texture->SetPixel(j, i, color);
}

static void CreateSpacecraftTexture() {
    spacecraft_emission_map = new sreTexture(256, 128);
    Fill(spacecraft_emission_map, 0, 0, 256, 128, 0);
    unsigned int yellow = 150 | (150 << 8) | 0xFF000000;
    unsigned int bright_yellow = 255 | (255 << 8) | 0xFF000000;
    unsigned int grey = 150 | (150 << 8) | (150 << 16) | 0xFF000000;
    unsigned int red = 150 | 0xFF000000;
    Fill(spacecraft_emission_map, 0, 0, 256, 8, grey);     // Top circle.
    Fill(spacecraft_emission_map, 0, 120, 256, 8, bright_yellow); // Bottom circle.
    for (int i = 0; i < 16; i++) {
        Fill(spacecraft_emission_map, i * 16, 61, 6, 6, yellow); // Colored windows on the sides.
        Fill(spacecraft_emission_map, i * 16 + 8, 61, 6, 6, red);
    }
    for (int j = 0; j < 26; j += 4)
        for (int i = 0; i < 64; i++) {
            if ((i & 7) == 7)
                continue;
            Fill(spacecraft_emission_map, i * 4, 24 + j, 2, 2, yellow); // Small windows.
            Fill(spacecraft_emission_map, i * 4, 103 - j, 2, 2, yellow);
        }
    spacecraft_emission_map->type = TEXTURE_TYPE_LINEAR;
    spacecraft_emission_map->UploadGL(false);
}

static void PlanetCreateScene(sreScene *scene, sreView *view) {
    // Disable shadow volumes should improve performance a little
    // (shadow volumes work, but show some artifacts). 16-bit vertices
    // indices can be used for more mesh models, and some GPU memory is
    // freed.
    sreSetShadowVolumeSupport(false);
    int physics_flag = 0;
    if (!physics)
        physics_flag = SRE_OBJECT_NO_PHYSICS;
    int hidden_flag = 0;
    int shadows_flag = SRE_OBJECT_CAST_SHADOWS;
    if (!show_spacecraft) {
        hidden_flag = SRE_OBJECT_HIDDEN;
        shadows_flag = 0;
    }
    Vector3D initial_ascend_vector = Vector3D(1.0, 0, 0);
//    Vector3D initial_ascend_vector = Vector3D(0, 0, 1.0);
    initial_ascend_vector.Normalize();
    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
        scene->SetAmbientColor(Color(0.03f, 0.03f, 0.03f)); //DEBUG
    else
        scene->SetAmbientColor(Color(0.2, 0.2, 0.2));
    sreModel *globe_model = sreCreateSphereModel(scene, 0);
    if (create_spacecraft) {
        // Add player sphere as scene object 0.
        scene->SetFlags(SRE_OBJECT_DYNAMIC_POSITION | shadows_flag |
            SRE_OBJECT_USE_TEXTURE | hidden_flag);
        scene->SetTexture(sreCreateStripesTexture(TEXTURE_TYPE_LINEAR,
            256, 256, 32, Color(0, 0.5f, 0.8f), Color(0.9f, 0.9f, 1.0f)));
        scene->SetDiffuseReflectionColor(Color(1.0f, 1.0f, 1.0f));
        scene->SetMass(3.0f);
        if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
            Point3D pos = Point3D(0, 0, 0) + (PLANETARY_RADIUS + 200.0f) * initial_ascend_vector;
            player_object = scene->AddObject(globe_model, pos.x, pos.y, pos.z, 0, 0, 0, 3.0);
        }
        else
            player_object = scene->AddObject(globe_model, 0, 0, 200.0f, 0, 0, 0, 3.0);
    }
    // Add player spacecraft as scene object 1.
    Point3D spacecraft_pos;
    sreModel *spacecraft_model;
    if (!create_spacecraft)
        goto skip_spacecraft;
    spacecraft_model = sreCreateEllipsoidModel(scene, 1.0, 0.3);
    scene->SetDiffuseReflectionColor(Color(0.1, 0.5, 0.1));
    scene->SetSpecularReflectionColor(Color(0.2, 0.2, 0.2));
    CreateSpacecraftTexture();
    scene->SetEmissionColor(Color(1.0, 1.0, 1.0));  
    scene->SetEmissionMap(spacecraft_emission_map);
    scene->SetFlags(SRE_OBJECT_DYNAMIC_POSITION | shadows_flag |
        SRE_OBJECT_USE_EMISSION_MAP | hidden_flag);
    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
        spacecraft_pos = Point3D(0.0f, 0.0f, 0.0f) + (PLANETARY_RADIUS + 300.0f) * initial_ascend_vector;
        spacecraft_object = scene->AddObject(spacecraft_model, spacecraft_pos.x, spacecraft_pos.y, spacecraft_pos.z, 0, 0, 0, 8.0);
    }
    else
        spacecraft_object = scene->AddObject(spacecraft_model, 0, 0, 300.0, 0, 0, 0, 8.0);
skip_spacecraft :
    scene->SetSpecularReflectionColor(Color(1.0, 1.0, 1.0));

    // Add sun sphere.
    scene->SetFlags(SRE_OBJECT_DYNAMIC_POSITION | SRE_OBJECT_NO_PHYSICS |
        SRE_OBJECT_EMISSION_ONLY | SRE_OBJECT_INFINITE_DISTANCE);
    scene->SetEmissionColor(Color(3.0, 3.0, 2.4));
    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
        sun_object = scene->AddObject(globe_model, 1000000.0, 0, 0, 0, 0, 0, 30000.0);
    else {
        Vector3D sun_N = Vector3D(0.6, 0.8, 0.5).Normalize();
        sun_object = scene->AddObject(globe_model, sun_N.x * 1000000.0, sun_N.y * 1000000.0, sun_N.z * 1000000.0, 0, 0, 0,
             30000.0);
    }
    scene->SetEmissionColor(Color(0, 0, 0));

    // Add terrain.
    // Create the terrain models.
    sreModel **mesh_model;
    int nu_models;
    CreatePlanetMesh(scene,  elevation_map_filename, elevation_map_texture,
        elevation_map_detail_scale_down, color_map_texture, 1, 1,
        mesh_flags, nu_models, mesh_model);

    // Add the terrain to the scene.
    scene->SetDiffuseReflectionColor(Color(1.0, 1.0, 1.0));
    scene->SetTexture(color_map_texture);
    // Because the Earth shader always uses emission and specular maps,
    // they have to be defined.
    int emission_flag = 0;
    if (night_map_filename != NULL)
        emission_flag = SRE_OBJECT_USE_EMISSION_MAP;
    int specularity_flag = 0;
    if (specularity_map_filename != NULL)
        specularity_flag = SRE_OBJECT_USE_SPECULARITY_MAP;
    int texture_flag = SRE_OBJECT_USE_TEXTURE;
    int multi_color_flag = 0;
    int earth_shader_flag = SRE_OBJECT_EARTH_SHADER;
    if (mesh_flags & MESH_FLAG_VERTEX_COLORS) {
        texture_flag = 0;
        multi_color_flag = SRE_OBJECT_MULTI_COLOR;
        // The Earth shader (which dynamically adjusts light levels) assumes the use of a
        // terrain color texture and not color attributes, so it can't be used with the
        // --vertex-colors option.
        earth_shader_flag = 0;
    }
    else if ((mesh_flags & MESH_FLAG_SPHERICAL_MODEL) && night_map_texture == NULL)
        earth_shader_flag = 0;
    // Emission color must be set to black for earth shader.
    scene->SetEmissionColor(Color(0.0f, 0.0f, 0.0f));
    if (night_only) {
        // Use only night texture, do not use earth shader.
        scene->SetEmissionMap(night_map_texture);
        scene->SetEmissionColor(Color(1.0f, 1.0f, 1.0f));
        scene->SetFlags(emission_flag | SRE_OBJECT_USE_TEXTURE |
            SRE_OBJECT_USE_SPECULARITY_MAP | SRE_OBJECT_CAST_SHADOWS |
            physics_flag);
    }
    else {
        if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
            scene->SetFlags(texture_flag | multi_color_flag | earth_shader_flag |
                emission_flag | specularity_flag | SRE_OBJECT_CAST_SHADOWS | physics_flag);
        else
            scene->SetFlags(texture_flag | multi_color_flag | specularity_flag |
                SRE_OBJECT_CAST_SHADOWS | physics_flag);
    }

    sreTexture *black_texture = sreCreateCheckerboardTexture(TEXTURE_TYPE_NORMAL,
        64, 64, 64, 64, Color(0, 0, 0), Color (0, 0, 0));
    if (specularity_map_filename != NULL) {
        scene->SetSpecularReflectionColor(Color(1.0f, 1.0f, 1.0f));
        scene->SetSpecularExponent(120.0f);
        scene->SetSpecularityMap(specularity_map_texture);
    }
    else {
        scene->SetSpecularReflectionColor(Color(0.0f, 0.0f, 0.0f));
        scene->SetSpecularExponent(1.0f);
        scene->SetSpecularityMap(black_texture);
    }

    if ((mesh_flags & MESH_FLAG_SPHERICAL_MODEL) && night_map_filename != NULL)
        // Set city light emission map for earth shader.
        scene->SetEmissionMap(night_map_texture);
    else
        scene->SetEmissionMap(black_texture);
    for (int i = 0; i < nu_models; i++) {
        // Dynamic LOD, starting from level 0 with threshold scale factor of 6.0;
        // there are one, two or three LOD levels, depending on the gains at
        // subsequent LOD levels; the second highest level is used for physics.
        int physics_lod_level = mesh_model[i]->nu_lod_levels - 2;
        physics_lod_level = maxi(0, physics_lod_level);
        scene->SetLevelOfDetail(SRE_LOD_DYNAMIC, 0, - 1,
            extra_lod_threshold_scaling, physics_lod_level);
        scene->AddObject(mesh_model[i], 0, 0, 0, 0, 0, 0, 1.0f);
    }
    scene->SetLevelOfDetail(SRE_LOD_DYNAMIC, 0, - 1, 1.0, 0);
    scene->SetSpecularReflectionColor(Color(1.0f, 1.0f, 1.0f));

    // Lights.
    // Add sun.
    if (!night_only) {
        Vector3D lightdir = Vector3D(- 0.6, - 0.8, - 0.5);
        lightdir.Normalize();
        directional_light = scene->AddDirectionalLight(SRE_LIGHT_DYNAMIC_DIRECTION,
            lightdir, Color(0.7f, 0.7f, 0.7f) * sun_light_factor);
    }
    if (night_only || (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)) {
        if (create_spacecraft) {
            Vector3D spot_dir = Vector3D(0, 1.0f, - 0.5f).Normalize();
            Vector3D rel_pos = Vector3D(0, 0, 0); // Vector3D(0, 0, - 2.4f);
            spacecraft_spot_light = scene->AddSpotLight(
                SRE_LIGHT_DYNAMIC_POSITION | SRE_LIGHT_DYNAMIC_DIRECTION,
                spacecraft_pos - rel_pos, spot_dir, 27.0f, 500.0f, Color(1.2f, 1.2f, 1.2f));
            scene->AttachLight(spacecraft_object, spacecraft_spot_light, Vector3D(0, 0, 0),
                spot_dir);
        }
    }

    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
        sre_internal_application->SetFlags(sre_internal_application->GetFlags() |
            SRE_APPLICATION_FLAG_DYNAMIC_GRAVITY | SRE_APPLICATION_FLAG_NO_GROUND_PLANE);
        sre_internal_application->gravity_position.Set(0, 0, 0);
    }
    else
        view->SetViewModeFollowObject(0, 40.0, Vector3D(0, 0, 10.0));
    if (create_spacecraft) {
        sre_internal_application->control_object = spacecraft_object;
        if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
            sre_internal_application->hovering_height = Magnitude(ProjectOnto(scene->object[
                sre_internal_application->control_object]->position,
                initial_ascend_vector));
        else
            sre_internal_application->hovering_height =
                scene->object[sre_internal_application->control_object]->position.z;
    }
    sreSetHDRKeyValue(0.2f);

    sre_internal_application->SetFlags(sre_internal_application->GetFlags() |
        SRE_APPLICATION_FLAG_NO_GRAVITY);
}

// #define SPECIFIC_PRECISION_DOUBLE

#ifdef SPECIFIC_PRECISION_DOUBLE
#define VectorSpecificPrecision3D VectorDouble3D
#define PointSpecificPrecision3D PointDouble3D
#else
#define VectorSpecificPrecision3D Vector3D
#define PointSpecificPrecision3D Point3D
#endif

static char message[80], message2[80];
static VectorDouble3D forward_vector;
static VectorDouble3D ascend_vector;
static VectorDouble3D right_vector;
#define HOUR_OFFSET 11.0f

void PlanetApplication::StepBeforeRender(double demo_time) {
    // Use double precision math to avoid unsmooth movement.
    PointSpecificPrecision3D sun_pos;
    if ((mesh_flags & MESH_FLAG_SPHERICAL_MODEL) && !night_only)
        goto skip_sun_pos;
    {
    double h = demo_time + HOUR_OFFSET * day_interval / 24.0d;
    double hour = fmod(h, day_interval) * 24.0d / day_interval;
    double day = floor(fmod(hour / 24.0d, 365.0d));

    if (seasonal_bmng) {
        int month = (int)floor(day * 12.0f / 365.0f);
        if (month != current_month) {
            current_month = month;
            seasonal_color_map_texture->opengl_id =
                bmng_seasonal_color_map_texture[month]->opengl_id;
        }
    }

    Matrix3D sr1;
    sr1.AssignRotationAlongZAxis(- hour * 2.0 * M_PI / 24.0);
    Matrix3D sr2;
    sr2.AssignRotationAlongXAxis(23.4 * M_PI / 180.0);
    Matrix3D sr3;
    sr3.AssignRotationAlongZAxis(fmod(demo_time, (day_interval * 365.0)) * 2.00 * M_PI /
        (day_interval * 365.0));
    sun_pos = PointSpecificPrecision3D(- PLANETARY_RADIUS * 1000.0d, 0, 0);
    sun_pos = PointSpecificPrecision3D(((sr3 * sr2) * sr1) * sun_pos.GetPoint3D());
    Vector3D light_dir = (- sun_pos).Normalize().GetVector3D();
    scene->ChangeDirectionalLightDirection(directional_light, light_dir);
    scene->ChangePosition(sun_object, sun_pos.x, sun_pos.y, sun_pos.z);
    if (display_time) {
         sprintf(message, "%02d:%02dh Day %d", (int)floorf(hour),
            (int)floorf((hour - (int)floorf(hour)) * 60.0 / 100.0),
            (int)(floor(day) + 1));
        sre_internal_application->text_message[0] = message; 
        sre_internal_application->text_message_time = sre_internal_backend->GetCurrentTime();
        sre_internal_application->nu_text_message_lines = 1;
    }
    }
skip_sun_pos :
    if (!create_spacecraft)
        return;

    saved_hovering_height = sre_internal_application->hovering_height;

    ascend_vector = VectorSpecificPrecision3D(0.0d, 0.0d, 1.0d);
    float view_distance;
    if (sre_internal_application->control_object == spacecraft_object) {
        if (show_spacecraft)
            view_distance = 100.0;
        else
            view_distance = 0.1f;
        // Hide the player (ball) object.
        scene->object[player_object]->flags |= SRE_OBJECT_HIDDEN;
        scene->object[player_object]->flags &= ~SRE_OBJECT_CAST_SHADOWS;
        // Hide the spacecraft object when required.
        if (!show_spacecraft) {
            scene->object[spacecraft_object]->flags |= SRE_OBJECT_HIDDEN;
            scene->object[spacecraft_object]->flags &= ~SRE_OBJECT_CAST_SHADOWS;
        }
    }
    else {
        view_distance = 40.0;
        // Show the player (ball) object.
        scene->object[player_object]->flags &= ~SRE_OBJECT_HIDDEN;
        scene->object[player_object]->flags |= SRE_OBJECT_CAST_SHADOWS;
        // Also show the spacecraft up in the air.
        scene->object[spacecraft_object]->flags &= ~SRE_OBJECT_HIDDEN;
        scene->object[spacecraft_object]->flags |= SRE_OBJECT_CAST_SHADOWS;
    }
    if (!(mesh_flags & MESH_FLAG_SPHERICAL_MODEL)) {
        view->SetViewModeFollowObject(sre_internal_application->control_object, view_distance,
            Vector3D(0, 0, 0) /* Vector3D(0, 0, view_distance * 0.25) */);
        return;
    }
    // Set viewing direction.
    VectorSpecificPrecision3D up_vector = VectorSpecificPrecision3D(scene->object[
        sre_internal_application->control_object]->position.GetVector3D());
    up_vector.Normalize();
    sre_internal_application->view->SetAscendVector(up_vector.GetVector3D());
    ascend_vector = up_vector;
    // Define the basal forward direction as looking down a meridian from the north pole (negative latitude direction).
    double latitude = asin(ascend_vector.z);
    double longitude = atan2(ascend_vector.y, ascend_vector.x);
    float height = Magnitude(scene->object[sre_internal_application->control_object]
        ->position.GetVector3D()) - PLANETARY_RADIUS;
    sprintf(message2, "%.2lf%c %.2lf%c H = %.0f", fabs(latitude * 180.0f / M_PI), latitude < 0 ? 'S' : 'N',
        fabs(longitude * 180.0f / M_PI), longitude < 0 ? 'W' : 'E', height);
    sre_internal_application->text_message[1] = message2;
    sre_internal_application->nu_text_message_lines = 2;

    Vector3D angles;
    sre_internal_application->view->GetViewAngles(angles);

#if 0
    // Rotate the forward direction according to the viewing angle.
    MatrixDouble3D r2; 
    r2.AssignRotationAlongAxis(ascend_vector, - angles.x * M_PI / 180.0);
    VectorDouble3D up_vector2 = VectorSpecificPrecision3D(r2 * up_vector);
#endif

#if 0
    // Calculate the great circle corresponding to the viewing direction.
    GreatCircleRoute great_circle;
    great_circle.lat0 = latitude;
    great_circle.long0 = longitude;
    great_circle.height0 = height;
    great_circle.height1 = height;
    great_circle.pitch0 = angles.x;
    great_circle.pitch1 = angles.x;

    // Calculate the view direction by moving a very small amount along the great circle to point P.
    double P_longitude, P_latitude, azimuth;
    great_circle.CalculateIntermediatePoint(0.001d, P_longitude, P_latitude, azimuth);
    // Calculate view vector.
    Matrix3D m1, m2;
    m1.AssignRotationAlongZAxis(P_longitude);
    m2.AssignRotationAlongYAxis(- P_latitude);
    Point3D target = m1 * (m2 * Point3D(1.0f, 0.0f, 0.0f));
    Vector3D view_direction = (target - viewpoint).Normalize();
    Vector3D right_vector = Cross(view_direction, up_vector).Normalize();
    Matrix3D m3;
    m3.AssignRotationAlongAxis(right_vector, - angles.x * M_PI / 180.0);
    view_direction = m3 * view_direction;

    Point3D viewpoint = (VectorSpecificPrecision3D(scene->object[
        sre_internal_application->control_object]->position) - (view_distance * view_direction)).
        GetVector3D();
    // + up_vector * view_distance * 0.25;

    sre_internal_application->view->SetForwardVector(forward_vector.GetVector3D());
    Point3D lookat = scene->object[sre_internal_application->control_object]->position;
    up_vector = VectorSpecificPrecision3D(m3 * up_vector);
#endif

#if 1
        // Define two arbitrary points on the great circle defined by the up_vector and thetaz (angles.z).
        // The base direction is negative latitude.
        double latitude1 = latitude - M_PI / 4.0d;
        double longitude1 = longitude;
        if (latitude1 <  - 0.5d * M_PI) {
            latitude1 = - M_PI - latitude1;
            longitude1 -= M_PI;
        }
        PointSpecificPrecision3D P1 = PointSpecificPrecision3D(cos(latitude1) * cos(longitude1), cos(latitude1) *
            sin(longitude1), sin(latitude1));
        double latitude2 = latitude + M_PI / 4.0d;
        double longitude2 = longitude;
        if (latitude2 > 0.5d * M_PI) {
            latitude2 = M_PI - latitude2;
            longitude2 += M_PI;
        }
        PointSpecificPrecision3D P2 = PointSpecificPrecision3D(cos(latitude2) * cos(longitude2), cos(latitude2) *
            sin(longitude2), sin(latitude2));
        MatrixDouble3D r1;
        double theta = angles.z * M_PI / 180.0d;
        r1.AssignRotationAlongAxis(up_vector, theta);
        P1 = PointSpecificPrecision3D(r1 * P1);
        P2 = PointSpecificPrecision3D(r1 * P2);
        // Calculate the normal of the new great circle.
        VectorSpecificPrecision3D great_circle_normal = VectorSpecificPrecision3D(CalculateNormal(
            PointSpecificPrecision3D(up_vector), P1, P2));
       
#ifdef DEBUG_MATH
        char *p1_str = P1.GetString();
        char *p2_str = P2.GetString();
        char *normal_str = great_circle_normal.GetString();
        sreMessage(SRE_MESSAGE_INFO, "%s%s%s%s%s%s", "P1 = ", p1_str, " P2 = ", p2_str,
            " Great circle normal = ", normal_str);
        delete [] p1_str;
        delete [] p2_str;
        delete normal_str;
#endif
//        if (longitude < 0)
//            great_circle_normal = - great_circle_normal;
//        if (latitude < 0)
//            great_circle_normal = - great_circle_normal;

    right_vector = great_circle_normal;
    forward_vector = Cross(right_vector, up_vector);
    forward_vector.Normalize();

#if 0
    Matrix3D r1;
    r1.AssignRotationAlongAxis(up_vector, angles.z * M_PI / 180);
    forward_vector = r1 * forward_vector;
    Vector3D right_vector = Cross(forward_vector, up_vector);
    right_vector.Normalize();
#endif
    sre_internal_application->view->SetForwardVector(forward_vector.GetVector3D());
    MatrixDouble3D r2; 
    r2.AssignRotationAlongAxis(right_vector.GetVector3D(), - angles.x * M_PI / 180.0);
    VectorSpecificPrecision3D view_direction = VectorSpecificPrecision3D(r2 * forward_vector.GetVector3D());
    Point3D viewpoint = (VectorSpecificPrecision3D(scene->object[
        sre_internal_application->control_object]->position) - (view_distance * view_direction)).
        GetVector3D();
    // + up_vector * view_distance * 0.25;
#if 0
    Point3D lookat = scene->object[sre_internal_application->control_object]->position; // viewpoint + view_direction;
#else
    Point3D lookat = viewpoint + view_direction;
#endif
    up_vector = VectorSpecificPrecision3D(r2 * up_vector);

#endif

    sre_internal_application->view->SetViewModeLookAt(viewpoint, lookat, up_vector.GetVector3D());
    sre_internal_application->view->SetMovementMode(SRE_MOVEMENT_MODE_USE_FORWARD_AND_ASCEND_VECTOR);
    // Let spacecraft spot light point the right way.
    // This is now handled automatically by libsre.
//    scene->ChangeSpotOrBeamLightDirection(spacecraft_spot_light, - ascend_vector);
}

void PlanetApplication::StepBeforePhysics(double demo_time) {
    if (sre_internal_application->flags & SRE_APPLICATION_FLAG_NO_GRAVITY) {
        if (sre_internal_application->control_object == player_object) {
            sre_internal_application->hovering_height = saved_hovering_height;
            sre_internal_application->hovering_height_acceleration = 0;
            scene->BulletChangeVelocity(spacecraft_object, Vector3D(0, 0, 0));
        }
        sre_internal_application->control_object = spacecraft_object;
    }
    else {
        if (sre_internal_application->control_object == spacecraft_object) {
            // Drop the player from the spacecraft.
            scene->BulletChangeVelocity(spacecraft_object, Vector3D(0, 0, 0));
            Point3D new_pos = (VectorDouble3D(scene->object[spacecraft_object]->position) -
                (ascend_vector * 15.0d)).GetVector3D();
//            scene->ChangePosition(player_object, new_pos);
            scene->BulletChangeVelocity(player_object, Vector3D(0, 0, 0));
            scene->BulletChangePosition(player_object, new_pos);
        }
        sre_internal_application->control_object = player_object;
    }
    MatrixDouble3D spin_matrix;
    if (show_spacecraft)
        spin_matrix.AssignRotationAlongZAxis(fmod(demo_time, 4.0d) * 2.0d * M_PI / 4.0d);
    else {
        spin_matrix.SetIdentity();
    }
    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
        // Try to keep the spacecraft upright, parallel to the surface.
        if (sre_internal_application->control_object == spacecraft_object) {
            MatrixDouble3D rot_matrix;
            rot_matrix.Set(ascend_vector.GetVector3D(), forward_vector.GetVector3D(), right_vector.GetVector3D());
            MatrixDouble3D r;
            r.AssignRotationAlongYAxis(M_PI / 2.0d);
//            rot_matrix.SetRow(0, right_vector);
//            rot_matrix.SetRow(1, forward_vector);
//            rot_matrix.SetRow(2, ascend_vector);
            saved_spacecraft_rotation_matrix = (rot_matrix * r) * spin_matrix;         
            scene->BulletChangeRotationMatrix(sre_internal_application->control_object,
               saved_spacecraft_rotation_matrix);
        }
        else
            // The spacecraft is not being controlled, but should rotate.
            scene->BulletChangeRotationMatrix(spacecraft_object, (saved_spacecraft_rotation_matrix *
                spin_matrix));
    }
    else
        scene->BulletChangeRotationMatrix(spacecraft_object, spin_matrix);
    // Set the maximum horizontal velocity (over the surface), for the spacecraft it increases as the height increases.
    if (sre_internal_application->control_object == player_object)
        sre_internal_application->max_horizontal_velocity = 100.0f;
    else {
        float height = Magnitude(ProjectOnto(scene->object[spacecraft_object]->position, ascend_vector.GetVector3D()));
        if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL)
            height -= PLANETARY_RADIUS;
        sre_internal_application->max_horizontal_velocity = 5.0f + height * 0.005f;
        sre_internal_application->horizontal_acceleration =
            sre_internal_application->max_horizontal_velocity * 2.0f;
        // The ascend/descend controls are also sensitive to the height above the surface.
        sre_internal_application->hovering_height_acceleration = 100.0f + height * 0.5f;
    }
    if (mesh_flags & MESH_FLAG_SPHERICAL_MODEL) {
        // Set viewing distance, clip distance increases as height increases.
        float player_dist = Magnitude(scene->object[sre_internal_application->control_object]->position) -
            PLANETARY_RADIUS;
        float far_plane_dist = maxf(2000.0, PLANETARY_RADIUS * 0.2f +
            player_dist * PLANETARY_RADIUS / 2500.0f);
        sreSetFarPlaneDistance(far_plane_dist);
        // Also increase the shadow mapping region as height increases.
        float factor;
        if (sre_internal_application->control_object == spacecraft_object)
            factor = far_plane_dist / 2000.0 + powf((far_plane_dist - 2000.0) / 2000.0, 0.2) * 4.0;
        else
            factor = 0.5;
        sreSetShadowMapRegion(Point3D(- 1000.0f, - 1000.0f, - 1000.0f) * factor,
            Point3D(1000.0f, 1000.0f, 200.0f) * factor);
//        sreSetShadowMapRegion(Point3D(- 400.0, - 400.0, - 600.0) * factor, Point3D(400.0, 400.0, 200.0) * factor);
    }
}

void PlanetSetParameters(float interval, bool _display_time, bool _physics,
bool _create_spacecraft, bool _show_spacecraft, float _sun_light_factor,
float _extra_lod_threshold_scaling) {
    day_interval = interval;
    display_time = _display_time;
    physics = _physics;
    create_spacecraft = _create_spacecraft;
    show_spacecraft = _show_spacecraft;
    sun_light_factor = _sun_light_factor;
    extra_lod_threshold_scaling = _extra_lod_threshold_scaling;
}

static void NoopTextOverlayFunc() {
}

enum {
    OPTION_ELEVATION_MAP = 0,
    OPTION_COLOR_MAP,
    OPTION_SPECULARITY_MAP,
    OPTION_NIGHT_MAP,
    OPTION_EARTH,
    OPTION_MOON,
    OPTION_MARS,
    OPTION_MERCURY,
    OPTION_EARTH_BMNG,
    OPTION_TEXTURE_SIZE,
    OPTION_ELEVATION_MAP_SCALE_DOWN,
    OPTION_NIGHT_ONLY,
    OPTION_NEW_CACHE,
    OPTION_NO_CACHE,
    OPTION_NO_TEXT_OVERLAY,
    OPTION_SEASONAL_BMNG,
    OPTION_VERTEX_COLORS,
    OPTION_NO_COLLISIONS,
};

static const char *option_string[][2] = {
    { "--elevation-map", "-e" },
    { "--color-map", "-c" },
    { "--specular-map", "-s" },
    { "--night-map", "-n" },
    { "--earth", "-E" },
    { "--moon", "-L" },
    { "--mars", "-M" },
    { "--mercury", "-Q" },
    { "--earth-bmng", "-bmng" },
    { "--texture-size", "-t" },
    { "--elevation-scale-down", "-d" },
    { "--night-only", "-o" },
    { "--update-cache", "-U" },
    { "--no-cache", "-C" },
    { "--no-text-overlay", "-T" },
    { "--seasonal-bmng", "-S" },
    { "--vertex-colors", "-v" },
    { "--no-collisions", "-O" },
};

static const char *option_description[] = {
    "Set elevation map texture.",
    "Set color map texture.",
    "Set specular map texture.",
    "Set night map texture.",
    "Set standard Earth texture set.",
    "Set Moon texture set.",
    "Set Mars texture set.",
    "Set Mercury texture set.",
    "Set BMNG Earth texture set.",
    "Limit texture size dimensions to n.",
    "Scale down elevation map (terrain mesh) detail by n (power of two), "
    "default 8 (1 is highest detail).",
    "Only display night texture (no sunlight).",
    "Update terrain model cache (force calculation).",
    "Do not use terrain model cache (always calculate).",
    "Do not display GUI text overlay (FPS, time, position etc).",
    "Use seasonal Earth BMNG data set textures (change every month).",
    "Store terrain vertex colors in GPU memory instead of terrain texture.",
    "No collision detection (Bullet).",
};

#define NU_OPTIONS (sizeof(option_string) / sizeof(option_string[0]))

void PlanetApplicationError(const char *format, ...) {
    va_list args;
    va_start(args, format);
    printf("(sre-planet) ");
    vprintf(format, args);
    va_end(args);
    printf("\n");
    fflush(stdout);
    sreFinalizeApplication(sre_internal_application);
    exit(1);
}

static char *interpret_filename(char *filename) {
    if (strcmp(filename, "none") == 0)
        return NULL;
    return filename;
}

static void PrintUsage() {
    sreMessage(SRE_MESSAGE_INFO,
       "sre-planet -- planet renderer using SRE rendering engine\n"
       "\nOptions:\n");
    for (int i = 0; i < NU_OPTIONS; i++) {
       sreMessage(SRE_MESSAGE_INFO, "%s, %s", option_string[i][0], option_string[i][1]);
       sreMessage(SRE_MESSAGE_INFO, "    %s", option_description[i]);
    }
}

int main(int argc, char **argv) {
    if (argc == 1) {
        PrintUsage();
        exit(0);
    }
    PlanetApplication *app = new PlanetApplication;
    sreInitializeApplication(app, &argc, &argv);
    scene = app->scene;

    if (argc == 1) {
        PrintUsage();
    }

    elevation_map_filename = NULL;
    color_map_filename = NULL;
    specularity_map_filename = NULL;
    night_map_filename =  NULL;
    texture_size = 0;
    elevation_map_detail_scale_down = 8;
    night_only = false;
    mesh_flags = MESH_FLAG_SPHERICAL_MODEL;
    no_text_overlay = false;
    seasonal_bmng = false;
    planet_preset = - 1;

    for (int argi = 1; argi < argc; argi++) {
        int option = - 1;
        for (int i = 0; i < NU_OPTIONS; i++)
            for (int j = 0; j < 2; j++)
                if (strcmp(argv[argi], option_string[i][j]) == 0)
                    option = i;
        if (option < 0)
            PlanetApplicationError("Invalid option %s.", argv[argi]);
        // Handle options without an argument.
        switch (option) {
        case OPTION_EARTH :
        case OPTION_MOON :
        case OPTION_MERCURY :
        case OPTION_EARTH_BMNG : 
            planet_preset = option - OPTION_EARTH;
            continue;
        case OPTION_MARS :
            PlanetApplicationError("Option %s not yet implemented.", argv[argi]);
            continue;
        case OPTION_NIGHT_ONLY :
            night_only = true;
            continue;
        case OPTION_NEW_CACHE :
            mesh_flags |= MESH_FLAG_NEW_CACHE;
            continue;
        case OPTION_NO_CACHE :
            mesh_flags |= MESH_FLAG_NO_CACHE;
            continue;
        case OPTION_NO_TEXT_OVERLAY :
            no_text_overlay = true;
            continue;
        case OPTION_SEASONAL_BMNG :
            seasonal_bmng = true;
            continue;
        case OPTION_VERTEX_COLORS :
            mesh_flags |= MESH_FLAG_VERTEX_COLORS;
            continue;
        case OPTION_NO_COLLISIONS :
            physics = false;
            continue;
        }
        // Handle options with an argument.
        argi++;
        if (argi >= argc)
            PlanetApplicationError("Expected argument for option %s.",
                argv[argi - 1]);
        switch (option) {
        case OPTION_ELEVATION_MAP :
            elevation_map_filename = interpret_filename(argv[argi]);
            continue;
        case OPTION_COLOR_MAP :
            color_map_filename = interpret_filename(argv[argi]);
            continue;
        case OPTION_SPECULARITY_MAP :
            specularity_map_filename = interpret_filename(argv[argi]);
            continue;
        case OPTION_NIGHT_MAP :
            night_map_filename = interpret_filename(argv[argi]);
            continue;
        case OPTION_TEXTURE_SIZE :
            texture_size = atoi(argv[argi]);
            int i;
            for (i = 7; i <= 30; i++)
                if (texture_size == (1 <<  i))
                    break;
            if (i > 30)
                PlanetApplicationError(
                    "Invalid value specified for option --texture-size "
                    "(must be power of two in range 1024 to 65536).");
            continue;
        case OPTION_ELEVATION_MAP_SCALE_DOWN :
            elevation_map_detail_scale_down = atoi(argv[argi]);
            if (elevation_map_detail_scale_down < 1)
                PlanetApplicationError("Invalid value specified for option "
                    "--elevation-scale-down.");
            continue;
        }
    }
    // When either elevation map or color map is not specified, use the default
    // planet.
    if ((elevation_map_filename == NULL || color_map_filename == NULL) &&
    planet_preset < 0)
        planet_preset = PLANET_EARTH;
    if (planet_preset >= 0) {
	if (elevation_map_filename == NULL)
            elevation_map_filename = planet_preset_filename[planet_preset][0];
        if (color_map_filename == NULL)
            color_map_filename = planet_preset_filename[planet_preset][1];
        if (specularity_map_filename == NULL)
            specularity_map_filename = planet_preset_filename[planet_preset][2];
        if (night_map_filename == NULL)
            night_map_filename = planet_preset_filename[planet_preset][3];
    }

    elevation_map_texture = new sreTexture(
        elevation_map_filename, TEXTURE_TYPE_NORMAL | SRE_TEXTURE_TYPE_FLAG_NO_UPLOAD);
    if (elevation_map_texture->width % elevation_map_detail_scale_down != 0
    || elevation_map_texture->height % elevation_map_detail_scale_down != 0)
        PlanetApplicationError("Elevation map detail division factor too large.");

    if (texture_size == 0) {
        if (mesh_flags & MESH_FLAG_VERTEX_COLORS)
            color_map_texture = sreCreateTexture(color_map_filename,
                SRE_TEXTURE_TYPE_FLAG_USE_UNCOMPRESSED_TEXTURE | SRE_TEXTURE_TYPE_FLAG_NO_UPLOAD);
        else
            color_map_texture = sreCreateTexture(color_map_filename, TEXTURE_TYPE_NORMAL);
        texture_size = color_map_texture->width;
    }
    else {
        if (seasonal_bmng) {
            for (int i = 0; i < 12; i++) {
                char filename[256];
                sprintf(filename, "BMNG/world.2004%02d-16384x8192", i);
                bmng_seasonal_color_map_texture[i] = sreCreateTextureLimitLevelWidth(filename,
                    TEXTURE_TYPE_NORMAL, texture_size);
            }
            seasonal_color_map_texture = new sreTexture;
            *seasonal_color_map_texture = *bmng_seasonal_color_map_texture[0];
        }
        else {
            if (mesh_flags & MESH_FLAG_VERTEX_COLORS)
                color_map_texture = sreCreateTexture(color_map_filename,
                    SRE_TEXTURE_TYPE_FLAG_USE_UNCOMPRESSED_TEXTURE | SRE_TEXTURE_TYPE_FLAG_NO_UPLOAD);
            else
                color_map_texture = sreCreateTextureLimitLevelWidth(color_map_filename,
                    TEXTURE_TYPE_NORMAL, texture_size);
            if (texture_size > color_map_texture->width)
                PlanetApplicationError("Texture size must be less or equal than color map texture.");
        }
    }
    if (specularity_map_filename != NULL) {
        specularity_map_texture = sreCreateTextureLimitLevelWidth(specularity_map_filename,
            TEXTURE_TYPE_SPECULARITY_MAP, texture_size);
        if (texture_size > specularity_map_texture->width)
           PlanetApplicationError("Texture size must be less or equal than specularity map texture.");
    }
    if (night_map_filename != NULL) {
        night_map_texture = sreCreateTextureLimitLevelWidth(
            night_map_filename, TEXTURE_TYPE_NORMAL, texture_size);
        if (texture_size > night_map_texture->width)
            PlanetApplicationError("Texture size must be less or equal than night map texture.");
    }
    if (no_text_overlay)
        sreSetDrawTextOverlayFunc(NoopTextOverlayFunc);

    sreMessage(SRE_MESSAGE_INFO, 
       "Elevation map size %dx%d, scale down factor %d.",
       elevation_map_texture->width, elevation_map_texture->height,
       elevation_map_detail_scale_down);
    sreMessage(SRE_MESSAGE_INFO,
       "Color map size %dx%d, actual size used %dx%d.",
       color_map_texture->width, color_map_texture->height, texture_size, texture_size / 2);

    PlanetCreateScene(scene, app->view);

    // Terrain color texture memory storage is no longer needed when colors are stored
    // as attributes.
    if (mesh_flags & MESH_FLAG_VERTEX_COLORS)
        delete color_map_texture;

    // Look east.
    app->view->SetViewAngles(Vector3D(- 30.0f, 0, 270.0f));
    current_month = - 1;

    app->SetFlags(app->GetFlags() & (~SRE_APPLICATION_FLAG_DISPLAY_FPS));
    sreRunApplication(app);

    sreFinalizeApplication(app);
    exit(0);
}
