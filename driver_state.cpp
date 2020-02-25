#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    //std::cout << "TODO: allocate and initialize state.image_depth."<< std::endl;
    //std::cout << "COMPLETED: allocation and initialization of state.image_color"
              << std::endl;

    // Allocate memory for image_color
    state.image_color = new pixel[height * width];
    // Initialize image_color
    for (unsigned int i = 0; i < (height * width); ++i) {
        state.image_color[i] = make_pixel(0,0,0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    int triangles;
    int vertex_index = 0;

    //allocate array of size 3 for data_geometry in rasterize_triangle function
    data_geometry* geoArray = new data_geometry[3];
    data_vertex vertexData;

    switch (type) {
        case render_type::triangle:
            // Find number of triangles for
            triangles = state.num_vertices / 3;

            for ( int i = 0; i < triangles; ++i) {
                for ( int j = 0; j < 3; ++j) {
                    geoArray[j].data = state.vertex_data + vertex_index;
                    vertex_index += state.floats_per_vertex;
                }
                for ( int k = 0; k < 3; ++k) {
                    vertexData.data = geoArray[k].data;
                    state.vertex_shader(vertexData, geoArray[k], state.uniform_data);
                }
                rasterize_triangle(state, (const data_geometry**)(&geoArray));
            }
            break;

        case render_type::indexed:
            break;

        case render_type::fan:
            break;

        case render_type::strip:
            break;

        default:
            std::cerr << "ERROR: invalid render_type!" << std::endl;
    }

    delete[] geoArray;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    int x[3], y[3];

    // Optimization to ignore pixels outside of image
    int min_x, min_y;
    int max_x, max_y;

    // variables for k0 k1 and k2 for barycentric coordinate calculation
    float k0[3], k1[3], k2[3];
    // variables for barycentric weight and total area of triangle
    float area;
    float bary_weight[3];

    // convert into homogeneous coordinates
    const float widthOverTwo = state.image_width / 2.0f;
    const float heightOverTwo = state.image_height / 2.0f;
    for (int i = 0; i < 3; ++i) {
        x[i] = (int)(widthOverTwo * (*in)[i].gl_Position[X] / (*in)[i].gl_Position[W]
                + (widthOverTwo - 0.5f));
        y[i] = (int)(heightOverTwo * (*in)[i].gl_Position[Y] / (*in)[i].gl_Position[W]
                + (heightOverTwo - 0.5f));
    }


    // Calculate barycentric weights for alpha beta and gamma

    // alpha
    k0[T_A] = x[T_B] * y[T_C] - x[T_C] * y[T_B];
    k1[T_A] = y[T_B] - y[T_C];
    k2[T_A] = x[T_C] - x[T_B];

    // beta
    k0[T_B] = x[T_C] * y[T_A] - x[T_A] * y[T_C];
    k1[T_B] = y[T_C] - y[T_A];
    k2[T_B] = x[T_A] - x[T_C];

    // gamma
    k0[T_C] = x[T_A] * y[T_B] - x[T_B] * y[T_A];
    k1[T_C] = y[T_A] - y[T_B];
    k2[T_C] = x[T_B] - x[T_A];

    // Find total area
    area = 0.5f * (k0[T_A] - (x[T_A] * y[T_C] - x[T_C] * y[T_A]) + k0[T_C]);

    // Find min/max coordinates for x & y
    min_x = state.image_width - 1;
    min_y = state.image_height - 1;
    max_x = 0;
    max_y = 0;
    for (int i = 0; i < 3; ++i) {
        min_x = std::min(min_x, x[i]);
        min_y = std::min(min_y, y[i]);

        max_x = std::max(max_x, x[i]);
        max_y = std::max(max_y, y[i]);

    }

    // Checking for out of bound values
    if (min_x < 0) { min_x = 0; }
    if (min_y < 0) { min_y = 0; }
    if (max_x > state.image_width) { max_x = state.image_width; }
    if (max_y > state.image_height) { max_y = state.image_height; }


    // Iterate through each pixel and calculate the barycentric weights for
    // each.
    for (int i = min_x; i < max_x + 1; i++) {
        for (int j = min_y; j < max_y + 1; j++) {
            for (int k = 0; k < 3; k++) {
                // Calculation is not done doing the iterative approach
                // We're multiplying every time to find the barycentric
                bary_weight[k] = .5f * (k0[k] + (k1[k] * i)
                    + (k2[k] * j)) / area;
            }

            if (inTriangle(bary_weight)) {
                // At some point this will need to be changed to get the
                // actual color of the pixel.
                state.image_color[i + j * state.image_width] =
                    make_pixel(255, 255, 255);
            }
        }
    }
}


bool inTriangle(float * bary_weight) {
    for (int i = 0; i < 3; i++) {
        if (bary_weight[i] < 0) {
            return false;
        }
    }

    return true;
}
