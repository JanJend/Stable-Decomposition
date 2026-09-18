#include <grlina/draw_hf.hpp>
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

struct Grid {
    std::vector<double> x_grid{-1,1}, y_grid{-2,2};
    std::vector<std::vector<int>> values;
};
}

// Inspect exactly the pixels passed to the PNG writer, without a PNG decoder.
int stbi_write_png(const char*, int w, int h, int components, const void* data, int stride) {
    require(components == 4 && stride == w*4);
    width = w; height = h;
    const auto* bytes = static_cast<const unsigned char*>(data);
    image.assign(bytes,bytes+stride*h);
    return 1;
}

int main() {
    try {
        using graded_linalg::save_hilbert_png;
        const std::array<unsigned char,3> white{255,255,255}, low{212,234,255}, high{0,0,0};
        Grid grid{{-1,1},{-2,2},{{0,1},{10,100}}};
        save_hilbert_png(grid,"unused.png",1,100);
        require(width == 625 && height == 545); // small grids leave room for labels
        require(colour(150,385) == white); // zero at bottom left; Y increases upwards
        require(colour(150,135) == low);
        require(colour(400,385) == graded_linalg::detail::hilbert_colour(0.5));
        require(colour(400,135) == high);

        grid.values = {{1,1},{1,1}};
        save_hilbert_png(grid,"unused.png",1,1);
        require(colour(150,135) == low);
        require(colour(500,135) == low && colour(500,385) == low); // constant colourbar
        save_hilbert_png(grid,"unused.png",1,100);
        require(colour(150,135) == low && colour(500,135) != low); // shared comparison scale

        grid.values = {{0,0},{0,0}};
        save_hilbert_png(grid,"unused.png",1,0);
        require(colour(150,135) == white && colour(500,135) == white);
        for (int y=85; y<445; ++y) for (int x=100; x<460; ++x)
            require(colour(x,y) == white);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
