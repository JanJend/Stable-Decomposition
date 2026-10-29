#include <grlina/invariants/hilbert_function.hpp>
#include <grlina/invariants/landscape.hpp>
#include <iostream>

namespace {
int width, height;
std::vector<unsigned char> image;

void require(bool condition) {
    if (!condition) throw std::runtime_error("Hilbert plot pixel check failed");
}

std::array<unsigned char,3> colour(int x, int y) {
    const auto i = (y*width+x)*4;
    return {image[i],image[i+1],image[i+2]};
}
}

// Inspect the actual pixels sent to the in-memory PNG encoder.
int stbi_write_png_to_func(stbi_write_func* write, void* context, int w, int h,
                          int components, const void* data, int stride) {
    require(components == 4 && stride == w*4);
    width = w; height = h;
    const auto* bytes = static_cast<const unsigned char*>(data);
    image.assign(bytes,bytes+stride*h);
    unsigned char encoded[]{137,'P','N','G'};
    write(context,encoded,4);
    return 1;
}

int main() {
    try {
        using namespace graded_linalg;
        const std::array<unsigned char,3> white{255,255,255}, low{212,234,255}, high{0,0,0};
        HilbertFunction<double,int> function{{-1,0,1},{-2,0,2},{{0,1,1},{10,100,100},{10,100,100}},100};
        const auto bytes = hilbert_function_image(function,400,300,1,100);
        require(bytes == std::vector<unsigned char>({137,'P','N','G'}));
        require(width == 665 && height == 485); // resolution counts heatmap pixels
        require(colour(150,335) == white); // Y increases upwards
        require(colour(150,135) == low);
        require(colour(400,335) == hilbert_colour(0.5));
        require(colour(400,135) == high);

        function.values = {{1,1,1},{1,1,1},{1,1,1}};
        function.maximum = 1;
        hilbert_function_image(function,400,300,1,1);
        require(colour(150,135) == low);
        require(colour(550,135) == low && colour(550,335) == low);
        hilbert_function_image(function,400,300,1,100);
        require(colour(150,135) == low && colour(550,135) != low);

        function.values = {{0,0,0},{0,0,0},{0,0,0}};
        function.maximum = 0;
        hilbert_function_image(function,400,300);
        require(colour(150,135) == white && colour(550,135) == white);
        for (int y=85; y<385; ++y) for (int x=100; x<500; ++x)
            require(colour(x,y) == white);

        // Unequal coordinate intervals must occupy unequal pixel widths.
        HilbertFunction<double,int> irregular{{0,1,10},{0,1,10},{{1,1,1},{100,100,100},{100,100,100}},100};
        hilbert_function_image(irregular,400,300,1,100);
        require(colour(120,200) == low && colour(200,200) == high);
        hilbert_function_image(irregular,vec<double>{0.5,2},vec<double>{0.5,2},400,300,1,100);
        require(colour(100,200) == low && colour(400,200) == high);

        // Explicit axes omit the birth at x=1; the image still samples the
        // full induced table, while cropping to exactly [0.5,2] x [0,2].
        Module<R2GradedSparseMatrix<int>> module(R2GradedSparseMatrix<int>(0,2,{}, {}, {{0,0},{1,0}}));
        module.hilbert_function_image({0.5,2},{0,2},400,300);
        require(colour(100,200) == low && colour(400,200) == high);
        // Landscape heights use a linear white-to-blue scale. Two free
        // summands distinguish layer 1 from layer 2 at the same coordinate.
        const Landscape landscape(Module<R2GradedSparseMatrix<int>>(
            R2GradedSparseMatrix<int>(0,2,{}, {}, {{0,0},{0.5,0.5}})), {0,0,1,1}, 0.25);
        const std::array<unsigned char,3> quarter{191,191,255}, blue{0,0,255};
        landscape.image(1,5,5);
        require(width == 270 && height == 190);
        require(colour(101,88) == quarter); // (0.25,0.25), height 0.25 of max 1
        require(colour(104,85) == blue); // (1,1), height 1
        landscape.image(2,5,5);
        require(colour(101,88) == white); // Second summand has not been born.
        require(colour(104,85) == blue); // Second-layer height is 0.5.
        landscape.image(3,5,5);
        for (int x=100; x<105; ++x) for (int y=85; y<90; ++y)
            require(colour(x,y) == white);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
