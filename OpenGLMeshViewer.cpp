#include <iostream>
#include <sstream>
#include <memory>
#include <boost/program_options.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <spdlog/spdlog.h>
#include "types.h"
#include "utils/utils.h"
#include "utils/glshader.h"
#include "utils/gldrawdata.h"
#include "utils/material.h"
#include "utils/light.h"
#include "sphere.h"
#include "trimesh.h"

using namespace std;
using namespace meshview;

namespace fs = boost::filesystem;

// application-wide vao
GLuint vao_;

// window dimensions
Vec2i window_size_g{ 1, 1 };

// the models and their material and xform
vector<TriMesh::Ptr> meshes_vector;
vector<Mat4r> mesh_xforms;
vector<Real> translations;
Material::Ptr sphere_material_g;
Real scale;

vector<Mat4r> mesh_xforms_reset;
Real camera_z_position{ 2 };
Real camera_reset{ 2 };

Real starting_xpos;
Real starting_ypos;

// scene lights
vector<Light::Ptr> lights_g;

void
ScaleObject(TriMesh::Ptr mesh, Mat4r& mat, Vec3r& center)
{
    Vec3r bmin{ kInfinity,kInfinity, kInfinity };
    Vec3r bmax = -bmin;
    mesh->GetBoundingBox(bmin, bmax);

    Real x_size = bmax[0] - bmin[0];
    Real y_size = bmax[1] - bmin[1];
    Real z_size = bmax[2] - bmin[2];
    Real ratio = 0;
    if (x_size > y_size && x_size > z_size)
    {
        ratio = scale / x_size;
    }
    else if (y_size > x_size && y_size > z_size)
    {
        ratio = scale / y_size;
    }
    else
    {
        ratio = scale / z_size;
    }
    mat.diagonal() = Vec4r{ ratio,ratio,ratio,1 };
    Mat4r scale_matrix{ Mat4r::Identity() };
    scale_matrix.diagonal() = Vec4r{ ratio,ratio,ratio,1 };
    Vec4r bmin4{ bmin[0],bmin[1],bmin[2],1 };
    Vec4r bmax4{ bmax[0], bmax[1], bmax[2], 1 };
    Vec4r scaled_bmin4 = scale_matrix * bmin4;
    Vec4r scaled_bmax4 = scale_matrix * bmax4;
    Vec3r scaled_bmin{ scaled_bmin4[0],scaled_bmin4[1],scaled_bmin4[2] };
    Vec3r scaled_bmax{ scaled_bmax4[0],scaled_bmax4[1],scaled_bmax4[2] };
    Vec3r diagonal = scaled_bmax - scaled_bmin;
    center = scaled_bmin + (0.5 * diagonal);
}

void
GetViewAndProjectionMatrices(glm::mat4& view_matrix, glm::mat4& proj_matrix)
{
    // compute aspect ratio
    float aspect = static_cast<float>(window_size_g[0]) /
        static_cast<float>(window_size_g[1]);
    view_matrix = glm::lookAt(glm::vec3(0, 0, camera_z_position), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
    proj_matrix = glm::perspective(static_cast<float>(60.0 * kDEGtoRAD), aspect,
        0.01f, 50.0f);
}

void
Display(Real glfw_time)
{
    // clear window
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // make sure we have a valid object
    if (!meshes_vector[0])
        return;

    // get view and projection matrices
    glm::mat4 view_matrix, proj_matrix;
    GetViewAndProjectionMatrices(view_matrix, proj_matrix);

    // fill GLDraw data for the meshes
    int current_index = 0;
    for (auto mesh_g : meshes_vector)
    {
        GLDrawData draw_data;
        draw_data.SetModelMatrix(EigenToGLM(mesh_xforms[current_index]));
        draw_data.SetViewMatrix(view_matrix);
        draw_data.SetProjectionMatrix(proj_matrix);
        draw_data.SetMaterial(sphere_material_g);
        draw_data.SetLights(lights_g);
        // draw the mesh
        mesh_g->DrawGL(draw_data);
        current_index++;
    }
}

void
WindowResizeCallback(GLFWwindow*/*window*/, int width, int height)
{
    // set viewport to occupy full canvas
    window_size_g = Vec2i{ width, height };
    glViewport(0, 0, width, height);
}

void
Reset_xforms()
{
    mesh_xforms = mesh_xforms_reset;
}

// brief Keyboard callback function, which is called everytime a key
// on the keyboard is pressed.
// param[in] window pointer to glfw window (unused)
// param[in] key key that was pressed/released/repeated etc.
// param[in] scancode scancode for the keyboard
// param[in] action one of GLFW_PRESS, GLFW_REPEAT or GLFW_RELEASE
// param[in] mods integer values specifying which modifier keys
//                 (Shift, Alt, Ctrl, etc.) are currently pressed
void
KeyboardCallback(GLFWwindow*/*window*/, int key, int /*scancode*/, int action,
    int /*mods*/)
{
    if (!meshes_vector[0])
        return;
    Real new_z_pos = camera_z_position * 0.9;
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        switch (key) {
        case GLFW_KEY_X:
            if (new_z_pos >= 0.01)
                camera_z_position = new_z_pos;
            else
                camera_z_position = 0.01;
            break;
        case GLFW_KEY_Z:
            camera_z_position *= 1.1;
            break;
        case GLFW_KEY_SPACE:
            camera_z_position = camera_reset;
            Reset_xforms();
            break;
        default:
            break;
        }
    }
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state != GLFW_PRESS)
    {
        return;
    }
    Real delta_xpos = starting_xpos - xpos;
    Real delta_ypos = starting_ypos - ypos;
    delta_xpos *= kPi;
    delta_ypos *= kPi;
    delta_xpos /= 180;
    delta_ypos /= 180;
    delta_xpos *= -0.2;
    delta_ypos *= -0.2;

    // Make x rotation matrix
    Mat4r x_rotation_matrix{ Mat4r::Identity() };
    Real c = cos(delta_xpos);
    Real s = sin(delta_xpos);
    x_rotation_matrix(0, 0) = c;
    x_rotation_matrix(2, 0) = -s;
    x_rotation_matrix(0, 2) = s;
    x_rotation_matrix(2, 2) = c;

    // Make y rotation matrix
    Mat4r y_rotation_matrix{ Mat4r::Identity() };
    Real cy = cos(delta_ypos);
    Real sy = sin(delta_ypos);
    y_rotation_matrix(1, 1) = cy;
    y_rotation_matrix(1, 2) = -sy;
    y_rotation_matrix(2, 1) = sy;
    y_rotation_matrix(2, 2) = cy;

    for (auto& trans : mesh_xforms)
    {
        trans = y_rotation_matrix * x_rotation_matrix * trans;
    }
    starting_xpos = xpos;
    starting_ypos = ypos;
}

void
mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)// && action == GLFW_PRESS)
    {
        glfwGetCursorPos(window, &starting_xpos, &starting_ypos);
    }
}

// Create the main OpenGL (GLFW) window
// param[in] width window width
// param[in] height window height
// param[in] window_name window name/title
// return pointer to created window
GLFWwindow*
CreateGLFWWindow(int width, int height, const std::string& window_name)
{
    // init glfw
    if (!glfwInit()) {
        spdlog::error("glfwInit failed");
        return nullptr;
    }

    // on the Codio box, we can only use 3.1 (GLSL 1.4). On macos, we
    // can use OpenGL 4.1 and only core profile
#if !defined(__APPLE__)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
#else
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // create main glfw window
    auto* window = glfwCreateWindow(width, height, window_name.c_str(),
        nullptr, nullptr);
    if (!window) {
        spdlog::error("glfwCreatewindow failed");
        return nullptr;
    }
    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        spdlog::error("glewInit failed");
        return nullptr;
    }

    // enable vsync
    glfwSwapInterval(1);

    // handle window resize events
    glfwSetFramebufferSizeCallback(window, WindowResizeCallback);

    // set keyboard callback function
    glfwSetKeyCallback(window, KeyboardCallback);

    // set the mouse button callback function
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // set the mouse cursor callback function
    glfwSetCursorPosCallback(window, cursor_position_callback);

    // get the current window size and call the resize callback to
    // update the viewport with the correct aspect ratio
    int window_width, window_height;
    glfwGetWindowSize(window, &window_width, &window_height);
    WindowResizeCallback(window, window_width, window_height);

    return window;
}


bool
ParseArguments(int argc, char** argv, std::vector<std::string>* mesh_names)
{
    namespace po = boost::program_options;
    po::options_description desc("options");
    try {
        desc.add_options()
            ("help,h", "print usage")
            ("mesh_name,m",
                po::value<vector<std::string>>(mesh_names)->multitoken(),
                "Mesh filenames");

        // parse arguments
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout << desc << endl;
            return false;
        }
        po::notify(vm);
    }
    catch (std::exception& e) {
        cout << desc << endl;
        spdlog::error("{}", e.what());
        return false;
    }
    catch (...) {
        cout << desc << endl;
        spdlog::error("Invalid arguments");
        return false;
    }
    return true;
}

//! \brief Main executable function
int
main(int argc, char** argv)
{
    std::vector<string> mesh_names;
    if (!ParseArguments(argc, argv, &mesh_names))
        return -1;

    // create the main GLFW window
    auto window = CreateGLFWWindow(1280, 720, "Olio - Meshes");
    if (!window)
        return -1;

    // create VAO
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // create phong material for the sphere
    Vec3r ambient{ 0, 0, 0 }, diffuse{ .8, .8, 0 }, specular{ .5, .5, .5 };
    Real shininess{ 50 };
    ambient = diffuse;
    auto material = std::make_shared<PhongMaterial>(ambient, diffuse, specular, shininess);

    auto glshader = make_shared<GLPhongShader>();
    if (!glshader->LoadShaders("../shaders/phong_vert.glsl",
        "../shaders/phong_frag.glsl")) {
        spdlog::error("Failed to load shaders.");
        return -1;
    }

    // set the gl shader for the material
    material->SetGLShader(glshader);

    uint num_meshes = mesh_names.size();
    translations.resize(num_meshes);
    scale = 2 / static_cast<Real>(num_meshes);
    if (num_meshes % 2 == 0)
    {
        // ceiling is the highest index in the lower half
        int ceiling = (static_cast<Real>(num_meshes) / 2) - 1;
        uint odds = 1;
        // do lower half of indices in reverse
        for (int i = ceiling; i >= 0; --i)
        {
            translations[i] = -static_cast<Real>(odds) / static_cast<Real>(num_meshes);
            odds += 2;
        }
        uint corresponding_lower_index = 0;
        // copy positives of lower half to upper half
        for (int i = static_cast<int>(num_meshes) - 1; i > ceiling; --i)
        {
            translations[i] = -translations[corresponding_lower_index];
            corresponding_lower_index++;
        }
    }
    else
    {
        uint middle = std::floor(num_meshes / 2);
        uint outer = num_meshes - 1;
        // Set first half up to middle
        for (uint i = 0; i < middle; ++i)
        {
            translations[i] = -(static_cast<Real>(outer - (i * 2)) / static_cast<Real>(num_meshes));
        }
        // Set middle to 0
        translations[middle] = 0;
        // Set second half
        for (uint i = 0; i < middle; ++i)
        {
            translations[outer - i] = -(translations[i]);
        }
    }
    size_t mesh_index = 0;
    // Load all models
    for (auto model : mesh_names)
    {
        // create a Sphere instance
        auto mesh_g = std::make_shared<TriMesh>();//std::make_shared<Sphere>();

        // create a Mat4r transform
        Mat4r xform{ Mat4r::Identity() };

        // load the mesh(es)
        fs::path filepath(model);

        if (!mesh_g->Load(filepath))
            spdlog::error("Failed to load mesh");

        // set sphere's material
        sphere_material_g = material;
        meshes_vector.push_back(mesh_g);
        Vec3r center;
        ScaleObject(mesh_g, xform, center);
        xform(0, 3) = translations[mesh_index];
        xform(1, 3) = 0 - center[1];
        xform(2, 3) = 0 - center[2];
        mesh_xforms.push_back(xform);
        mesh_index++;
        mesh_xforms_reset.push_back(xform);
    }
    // add point light 1
    auto point_light1 = make_shared<PointLight>(Vec3r{ 2, 2, 4 }, Vec3r{ 10, 10, 10 },
        Vec3r{ 0.01f, 0.01f, 0.01f });
    lights_g.push_back(point_light1);

    // add point light 2
    auto point_light2 = make_shared<PointLight>(Vec3r{ -1, -4, 1 }, Vec3r{ 7, 2, 2 },
        Vec3r{ 0.01f, 0.01f, 0.01f });
    lights_g.push_back(point_light2);

    // add point light 3
    auto point_light3 = make_shared<PointLight>(Vec3r{ -2, 4, 1 }, Vec3r{ 0, 5, 2 },
        Vec3r{ 0.01f, 0.01f, 0.01f });
    lights_g.push_back(point_light3);

    // main draw loop
    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            break;
        Display(glfwGetTime());
        glfwSwapBuffers(window);
        glfwPollEvents();
        glfwWaitEvents();
    }

    // clean up stuff
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
