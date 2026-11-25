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
        const std::array<unsigned char,3> white{255,255,255}, low = hilbert_colour(0), high = hilbert_colour(1);
        HilbertFunction<double,int> function{{-1,0,1},{-2,0,2},{{0,1,1},{10,100,100},{10,100,100}},100};
        const auto bytes = hilbert_function_image(function,400,300,1,100);
        require(bytes == std::vector<unsigned char>({137,'P','N','G'}));
        require(width == 665 && height == 415); // resolution counts heatmap pixels
        require(colour(150,295) == white); // Y increases upwards
        require(colour(150,95) == low);
        require(colour(400,295) == hilbert_colour(0.5));
        require(colour(400,95) == high);

        auto has_ink = [&](int left, int top, int w, int h) {
            for (int x=left; x<left+w; ++x) for (int y=top; y<top+h; ++y)
                if (colour(x,y) != white) return true;
            return false;
        };
        require(!has_ink(100,0,400,44)); // No automatic heading or maximum subtitle.
        require(!has_ink(100,385,400,30)); // No default X caption or footer legend.
        require(!has_ink(10,45,40,300)); // No default Y caption.
        hilbert_function_image(function,400,300,1,100,"","Radius","Codensity");
        require(has_ink(100,385,400,25));
        require(has_ink(10,45,40,300)); // Long Y captions are drawn vertically.
        require(colour(150,295) == white && colour(400,95) == high);

        function.values = {{1,1,1},{1,1,1},{1,1,1}};
        function.maximum = 1;
        hilbert_function_image(function,400,300,1,1);
        require(colour(150,95) == low);
        require(colour(550,95) == low && colour(550,295) == low);
        hilbert_function_image(function,400,300,1,100);
        require(colour(150,95) == low && colour(550,95) != low);

        function.values = {{0,0,0},{0,0,0},{0,0,0}};
        function.maximum = 0;
        hilbert_function_image(function,400,300);
        require(colour(150,95) == white && colour(550,95) == white);
        for (int y=45; y<345; ++y) for (int x=100; x<500; ++x)
            require(colour(x,y) == white);

        // Unequal coordinate intervals must occupy unequal pixel widths.
        HilbertFunction<double,int> irregular{{0,1,10},{0,1,10},{{1,1,1},{100,100,100},{100,100,100}},100};
        hilbert_function_image(irregular,400,300,1,100);
        require(colour(120,160) == low && colour(200,160) == high);
        hilbert_function_image(irregular,vec<double>{0.5,2},vec<double>{0.5,2},400,300,1,100);
        require(colour(100,160) == low && colour(400,160) == high);

        // Explicit axes omit the birth at x=1; the image still samples the
        // full induced table, while cropping to exactly [0.5,2] x [0,2].
        Module<R2GradedSparseMatrix<int>> module(R2GradedSparseMatrix<int>(0,2,{}, {}, {{0,0},{1,0}}));
        module.hilbert_function_image({0.5,2},{0,2},400,300,"Radius","Codensity");
        require(colour(100,160) == low && colour(400,160) == high);
        require(has_ink(100,385,400,25) && has_ink(10,45,40,300));
        // Landscape heights use the Hilbert palette on a linear scale. Two free
        // summands distinguish layer 1 from layer 2 at the same coordinate.
        const Landscape landscape(Module<R2GradedSparseMatrix<int>>(
            R2GradedSparseMatrix<int>(0,2,{}, {}, {{0,0},{0.5,0.5}})), {0,0,1,1}, 0.25);
        const auto quarter = hilbert_colour(0.25);
        landscape.image(1,5,5,ColourMap::blue,0,"r","t");
        require(width == 270 && height == 120);
        require(colour(101,48) == quarter); // (0.25,0.25), height 0.25 of max 1
        require(colour(104,45) == high); // (1,1), height 1
        landscape.image(2,5,5);
        require(colour(101,48) == white); // Second summand has not been born.
        require(colour(104,45) == high); // Second-layer height is 0.5.
        landscape.image(3,5,5);
        for (int x=100; x<105; ++x) for (int y=45; y<50; ++y)
            require(colour(x,y) == white);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
