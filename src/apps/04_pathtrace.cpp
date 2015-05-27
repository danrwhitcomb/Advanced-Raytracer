#include "scene.h"
#include "intersect.h"
#include "montecarlo.h"
#include "animation.h"

#include <thread>

using std::thread;

// modify the following line to disable/enable parallel execution of the pathtracer
bool parallel_pathtrace = true;

image3f pathtrace(Scene* scene, bool multithread);
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose);
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth);



// lookup texture value
vec3f lookup_scaled_texture(vec3f value, image3f* texture, vec2f uv, bool tile = false, bool bilinear_filter = false) {
    if(not texture) return value;
    
    float u, v;
    float u_prime, v_prime;
    int i, j;
    int width, height;
    vec3f color_1, color_2, color_3, color_4, filtered_color;

    width = texture->width();
    height = texture->height();

    if (tile) {
        // Tile the texture if texcoords exceed [0, 1)
        u = uv.x - floor(uv.x);
        v = uv.y - floor(uv.y);
    } else {
        // Simple clamp coordinates
        u = clamp(uv.x, 0.0f, 1.0f);
        v = clamp(uv.y, 0.0f, 1.0f);
    }

    // Return the simple value from one pixel if not bilinear filtering
    if (!bilinear_filter) {return value * texture->at(u*(width-1), v*(height-1));}

    // Bilinear filter
    i = floor(u * width);
    j = floor(v * height);

    // Retreive the pixels that we will interpolate with
    color_1 = texture->at(i-1, j-1);
    color_2 = texture->at(i, j-1);
    color_3 = texture->at(i-1, j);
    color_4 = texture->at(i, j);

    // Compute interpolation factors
    u_prime = width * u - floor(width * u);
    v_prime = height * v - floor(height * v);

    // Find the interpolated color vector
    filtered_color = (
                (1 - u_prime) * (1 - v_prime) * color_1 +
                (u_prime) * (1 - v_prime) * color_2 +
                (1 - u_prime) * (v_prime) * color_3 +
                u_prime * v_prime * color_4
    );

    return value * filtered_color;
}

// compute the brdf
vec3f eval_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f l, vec3f norm, bool microfacet) {
    if (not microfacet) {
        auto h = normalize(v+l);
        return kd/pif + ks*(n+8)/(8*pif) * pow(max(0.0f,dot(norm,h)),n);
    } else {
        put_your_code_here("Implement microfacet brdf");
        return one3f; // <- placeholder
    }
}

// evaluate the environment map
vec3f eval_env(vec3f ke, image3f* ke_txt, vec3f dir) {
    return ke; // <- placeholder
}


vec3f compute_blur_reflection(Scene* scene, vec3f kr, ray3f ray, int depth, Rng* rng, float bsz, int bsa){
    vec3f r_i = zero3f;
    vec3f u = vec3f(ray.d.y * -1.0, ray.d.x, 0);
    vec3f v = vec3f(0, -1.0 * ray.d.z, ray.d.y);
    
    for(int i=0; i < bsa; i++){
        vec2f rand = rng->next_vec2f();
        
        vec3f n = (ray.d +
                    (0.5f-rand.x)*bsz*u +
                    (0.5f-rand.y)*bsz*v);
        
        ray3f new_ray;
        new_ray.d = n / (length(n) * 1.0);
        new_ray.e = ray.e;
        
        r_i += pathtrace_ray(scene, new_ray, rng, depth+1);
    }
    
    return (kr * r_i) / (bsa * 1.0);
}

//Get cosine of angle between two vectors
float vector_cosine(vec3f a, vec3f b) { return dot(a, b) / (length(a) * length(b));}


// compute the color corresponing to a ray by pathtrace
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth) {

    // get scene intersection
    auto intersection = intersect(scene,ray);

    bool mipmap = intersection.mat->mipmap;
    
    // if not hit, return background (looking up the texture by converting the ray direction to latlong around y)
    if(not intersection.hit) {
        return eval_env(scene->background, scene->background_txt, ray.d);
    }
    
    
    // setup variables for shorter code
    auto pos = intersection.pos;
    auto norm = intersection.norm;
    auto v = -ray.d;

    vec3f ke, kd, ks;
    vec3f kd_1, kd_2, kd_3;

    bool tex_tile = intersection.mat->tex_tile;
    bool tex_filter = intersection.mat->tex_filter;
    
    // compute material values by looking up textures
    if (mipmap) {

        // Find how far the intersection was from the camera and limit to [1, 5]
        float dist_heuristic = dist(intersection.pos, scene->camera->frame.o);
        dist_heuristic = clamp(dist_heuristic, 1.0, 5.0);

        // Intersection is close: use just the high-res textures
        if (dist_heuristic == 1) {
            kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord, tex_tile, tex_filter);

        // Blend the high-res and medium-res textures
        } else if (dist_heuristic < 3) {
            float blend = (dist_heuristic - 1 ) / 2;
            kd_1 = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord, tex_tile, tex_filter);
            kd_2 = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt_2, intersection.texcoord, tex_tile, tex_filter);
            kd = (1 - blend) * kd_1 + blend * kd_2;

        // Blend the medium-res and low-res textures
        } else if (dist_heuristic < 5) {
            float blend = (dist_heuristic - 3) / 2;
            kd_2 = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt_2, intersection.texcoord, tex_tile, tex_filter);
            kd_3 = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt_3, intersection.texcoord, tex_tile, tex_filter);
            kd = (1 - blend) * kd_2 + blend * kd_3;

        // Intersection is far away: just use the low-res texture
        } else {
            kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt_3, intersection.texcoord, tex_tile, tex_filter);
        }
    }

    // Not using mipmap - use the base kd_txt to calculate diffuse
    else {
        kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord, tex_tile, tex_filter);
    }

    ke = lookup_scaled_texture(intersection.mat->ke, intersection.mat->ke_txt, intersection.texcoord, tex_tile, tex_filter);
    ks = lookup_scaled_texture(intersection.mat->ks, intersection.mat->ks_txt, intersection.texcoord, tex_tile, tex_filter);
    auto n = intersection.mat->n;
    auto mf = intersection.mat->microfacet;
    
    // accumulate color starting with ambient
    auto c = scene->ambient * kd;
    
    // add emission if on the first bounce
    if(depth == 0 and dot(v,norm) > 0) c += ke;
    
    // foreach point light
    for(auto light : scene->lights) {
        // compute light response
        auto cl = light->intensity / (lengthSqr(light->frame.o - pos));
        // compute light direction
        auto l = normalize(light->frame.o - pos);
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
        // multiply brdf and light
        auto shade = cl * brdfcos;
        // check for shadows and accumulate if needed
        if(shade == zero3f) continue;
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos,light->frame.o))) c += shade;
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // foreach surface
    for(Surface* surface : scene->surfaces){
        // skip if no emission from surface
        if(surface->mat->ke == zero3f){continue;}
        // todo: pick a point on the surface, grabbing normal, area, and texcoord
        // generate a 2d random number
        vec2f rand = rng->next_vec2f();
        vec3f pos, normal;
        float area;
        auto r = surface->radius;

        // check if quad
        if(surface->isquad){
            // compute light position, normal, area
            pos.x = ((r + r) * rand.x) - r;
            pos.y = ((r + r) * rand.y) - r;
            pos.z = 0;
            normal = surface->frame.z;
            area = pow(r+r, 2);
            
        }
        // else if sphere
        else {
            // compute light position, normal, area
            pos = sample_direction_spherical_uniform(rand);
            auto l = length(pos);
            pos = pos * (r/l);
            normal = transform_direction_from_local(surface->frame,pos);
            area = 4 * pi * pow(r, 2);
        }
        
        pos = transform_point_from_local(surface->frame, pos);
        
        // set tex coords as random value got before
        vec2f texcoord = rand;
        // get light emission from material and texture
        vec3f ke = lookup_scaled_texture(surface->mat->ke, surface->mat->ke_txt, texcoord);
        // compute light direction
        vec3f light_dir = normalize(pos - intersection.pos);
        // compute light response (ke * area * cos_of_light / dist^2)
        vec3f response = ke * area * vector_cosine(normal, -light_dir) / distSqr(intersection.pos, pos);
        // compute the material response (brdf*cos)
        auto mat = intersection.mat;
        auto brdfcos = max(dot(intersection.norm, light_dir),0.0f) * eval_brdf(mat->kd, mat->ks, mat->n, v, light_dir, intersection.norm, mf);
        // multiply brdf and light
        auto shade = brdfcos * response;
        
        // check for shadows and accumulate if needed
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos, intersection.pos))){
                c += shade;
            }
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // todo: sample the brdf for environment illumination if the environment is there
    // if scene->background is not zero3f
//    if(scene->background != zero3f && !intersection.hit){
//        // pick direction and pdf;
//        auto rand = rng->next_vec2f();
//        
//        
//        pair<vec3f, float> sample =
//        
//        // compute the material response (brdf*cos)
//        // todo: accumulate response scaled by brdf*cos/pdf
//        // if material response not zero3f
//            // if shadows are enabled
//                // perform a shadow check and accumulate
//                // else just accumulate
//    }
    
    // todo: sample the brdf for indirect illumination
    // if kd and ks are not zero3f and haven't reach max_depth
        // pick direction and pdf
        // compute the material response (brdf*cos)
        // accumulate recersively scaled by brdf*cos/pdf
    
    // if the material has reflections
    if(not (intersection.mat->kr == zero3f)) {
        auto rr = ray3f(intersection.pos,reflect(ray.d,intersection.norm));
        
        if(intersection.mat->bsz != 0 && intersection.mat->bsa != 0){
            c += compute_blur_reflection(scene, intersection.mat->kr, rr, depth, rng, intersection.mat->bsz, intersection.mat->bsa);
        } else {
            c += pathtrace_ray(scene, rr, rng, depth++);
        }
    }
    
    // return the accumulated color
    return c;
}


// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "04_pathtrace", "raytrace a scene",
            {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue() } },
            {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json") },
               {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("") } }
        });
    
    auto scene_filename = args.object_element("scene_filename").as_string();
    Scene* scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");
    
    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    
    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    
    // NOTE: acceleration structure does not support animations
    message("reseting animation...\n");
    animate_reset(scene);
    
    message("accelerating...\n");
    accelerate(scene);
    
    message("rendering %s...\n", scene_filename.c_str());
    auto image = pathtrace(scene, parallel_pathtrace);
    
    message("saving %s...\n", image_filename.c_str());
    write_png(image_filename, image, true);
    
    delete scene;
    message("done\n");
}


/////////////////////////////////////////////////////////////////////
// Rendering Code


// pathtrace an image
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose) {
    if(verbose) message("\n  rendering started        ");
    // foreach pixel
    for(auto j = offset_row; j < scene->image_height; j += skip_row ) {
        if(verbose) message("\r  rendering %03d/%03d        ", j, scene->image_height);
        for(auto i = 0; i < scene->image_width; i ++) {
            // init accumulated color
            image->at(i,j) = zero3f;
            // grab proper random number generator
            auto rng = &rngs->at(i, j);
            // foreach sample
            for(auto jj : range(scene->image_samples)) {
                for(auto ii : range(scene->image_samples)) {
                    // compute ray-camera parameters (u,v) for the pixel and the sample
                    auto u = (i + (ii + rng->next_float())/scene->image_samples) /
                        scene->image_width;
                    auto v = (j + (jj + rng->next_float())/scene->image_samples) /
                        scene->image_height;
                    // compute camera ray
                    auto ray = transform_ray(scene->camera->frame,
                        ray3f(zero3f,normalize(vec3f((u-0.5f)*scene->camera->width,
                                                     (v-0.5f)*scene->camera->height,-1))));
                    // set pixel to the color raytraced with the ray
                    image->at(i,j) += pathtrace_ray(scene,ray,rng,0);
                }
            }
            // scale by the number of samples
            image->at(i,j) /= (scene->image_samples*scene->image_samples);
        }
    }
    if(verbose) message("\r  rendering done        \n");
    
}

// pathtrace an image with multithreading if necessary
image3f pathtrace(Scene* scene, bool multithread) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);
    
    // create a random number generator for each pixel
    auto rngs = RngImage(scene->image_width, scene->image_height);

    // if multitreaded
    if(multithread) {
        // get pointers
        auto image_ptr = &image;
        auto rngs_ptr = &rngs;
        // allocate threads and pathtrace in blocks
        auto threads = vector<thread>();
        auto nthreads = thread::hardware_concurrency();
        for(auto tid : range(nthreads)) threads.push_back(thread([=](){
            return pathtrace(scene,image_ptr,rngs_ptr,tid,nthreads,tid==0);}));
        for(auto& thread : threads) thread.join();
    } else {
        // pathtrace all rows
        pathtrace(scene, &image, &rngs, 0, 1, true);
    }
    
    // done
    return image;
}


