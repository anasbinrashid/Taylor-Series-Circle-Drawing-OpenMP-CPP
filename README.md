# OpenMP Circle Rendering Performance Comparison

A parallel computing project that demonstrates the performance differences between serial and parallel implementations of trigonometric functions using OpenMP. The program compares standard math library functions with custom Taylor series implementations, both in serial and parallel versions.

## Project Overview

This assignment (PDC Assignment 03) creates four different circle renderings using:
- **Serial Math Library**: Standard `cos()` and `sin()` functions
- **Parallel Math Library**: OpenMP parallelized standard functions  
- **Serial Taylor Series**: Custom Taylor series implementation
- **Parallel Taylor Series**: OpenMP parallelized Taylor series

Each circle is rendered with different radii and colors to visualize the performance comparison results.

## Features

- **Performance Benchmarking**: Measures execution time for each implementation
- **Visual Rendering**: SDL2-based graphical output showing all four circles
- **Color-Coded Results**: Different colors for easy identification
- **OpenMP Parallelization**: Uses 4 threads with static scheduling
- **Taylor Series Implementation**: Custom mathematical approximations for trigonometric functions

## 🛠Dependencies

### Required Libraries
- **OpenMP**: For parallel processing
- **SDL2**: For graphics rendering and window management

### Installation Commands

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install libomp-dev libsdl2-dev
```

**macOS (using Homebrew):**
```bash
brew install libomp sdl2
```

**Windows (using vcpkg):**
```bash
vcpkg install openmp sdl2
```

## 🚀 Compilation and Execution

### Compile the Program
```bash
g++ -fopenmp -o openmp_circle OpenMP.cpp -lSDL2 -lm
```

### Run the Program
```bash
./openmp_circle
```

## Output Explanation

### Console Output
The program displays execution times for each implementation:
```
Serial Math Library Time: X.XXXX s
Parallel Math Library Time: X.XXXX s
Serial Taylor Series Time: X.XXXX s
Parallel Taylor Time: X.XXXX s
```

### Visual Output
A 800x800 pixel window displays four concentric circles:
- **White Circle**: Serial Math Library (radius: 200px)
- **Red Circle**: Parallel Math Library (radius: 195px)
- **Green Circle**: Serial Taylor Series (radius: 190px)
- **Blue Circle**: Parallel Taylor Series (radius: 185px)

## Configuration Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `screenwidth` | 800 | Window width in pixels |
| `screenheight` | 800 | Window height in pixels |
| `numberofpoints` | 1000 | Points per circle |
| `numberofterms` | 1500 | Taylor series terms |
| `omp_num_threads` | 4 | OpenMP thread count |

## Mathematical Implementation

### Taylor Series Formulations

**Cosine Series:**
```
cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + ...
```

**Sine Series:**
```
sin(x) = x - x³/3! + x⁵/5! - x⁷/7! + ...
```

### Parallelization Strategy
- Uses `#pragma omp parallel for` with static scheduling
- Applies parallelization to both Taylor series computation and circle point generation
- Fixed thread count of 4 for consistent performance measurement

## Expected Performance Results

Typical performance characteristics:
1. **Serial Math Library**: Fastest (optimized system functions)
2. **Parallel Math Library**: Minimal speedup (overhead > computation)
3. **Serial Taylor Series**: Slowest (complex computations)
4. **Parallel Taylor Series**: Significant speedup over serial Taylor

## Troubleshooting

### Common Issues

**Compilation Errors:**
- Ensure OpenMP is properly installed: `gcc -fopenmp --version`
- Check SDL2 installation: `pkg-config --libs sdl2`

**Runtime Issues:**
- **Black Screen**: Verify SDL2 libraries are properly linked
- **No Performance Gain**: Check if system supports multi-threading
- **Segmentation Fault**: Ensure sufficient memory for point arrays

### Performance Optimization Tips
- Adjust `numberofterms` for Taylor series accuracy vs. speed trade-off
- Modify `numberofpoints` for rendering quality vs. performance
- Experiment with different OpenMP scheduling policies (`dynamic`, `guided`)

## Learning Objectives

This assignment demonstrates:
- **OpenMP Fundamentals**: Parallel loop constructs and thread management
- **Performance Analysis**: Comparing serial vs. parallel execution times
- **Mathematical Computing**: Taylor series implementation and optimization
- **Graphics Programming**: SDL2 integration for visual output
- **Parallel Algorithm Design**: Understanding when parallelization provides benefits

## Contributing

Feel free to fork this repository and experiment with:
- Different parallelization strategies
- Alternative mathematical series implementations
- Enhanced visualization features
- Performance profiling tools
