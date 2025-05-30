#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <SDL2/SDL.h>

using namespace std;

const int screenwidth = 800;
const int screenheight = 800;
const int numberofpoints = 1000;
const int numberofterms = 1500;

double taylorseriescosine(double x) 
{
    double sum = 1.0;
    double term = 1.0;
    double xsquare = x*x;
    
    for (int i = 1; i <= numberofterms; ++i) 
    {
        term *= (-1) * xsquare / ((2 * i - 1) * (2 * i));
        sum += term;
    }
    
    return sum;
}

double paralleltaylorseriescosine(double x) 
{
    double sum = 1.0;
    double term = 1.0;
    double xsquare = x*x;
    omp_set_num_threads(4);
    #pragma omp parallel for schedule (static)
    for (int i = 1; i <= numberofterms; ++i) 
    {
        term *= (-1) * xsquare / ((2 * i - 1) * (2 * i));
        sum += term;
    }
    
    return sum;
}


double taylorseriessine(double x) 
{
    double sum = x;
    double term = x;
    double xsquare = x*x;
    
    for (int i = 1; i <= numberofterms; ++i) 
    {
        term *= (-1) * xsquare / ((2 * i) * (2 * i + 1));
        sum += term;
    }
    
    return sum;
}

double paralleltaylorseriessine(double x) 
{
    double sum = x;
    double term = x;
    double xsquare = x*x;
    omp_set_num_threads(4);
    #pragma omp parallel for schedule (static)
    for (int i = 1; i <= numberofterms; ++i) 
    {
        term *= (-1) * xsquare / ((2 * i) * (2 * i + 1));
        sum += term;
    }
    
    return sum;
}


int main() 
{
    double starttime, endtime;
    
    vector<SDL_Point> mathlibrary(numberofpoints);
    vector<SDL_Point> parallelmathlibrary(numberofpoints);

    starttime = omp_get_wtime();
    
    for (int t = 0; t < numberofpoints; t++) 
    {
        double radians = 2.0 * t * M_PI / numberofpoints;
        
        mathlibrary[t] = {static_cast<int>(200 * cos(radians) + screenwidth / 2), static_cast<int>(200 * sin(radians) + screenheight / 2)};
    }
    
    endtime = omp_get_wtime();
    
    cout << "Serial Math Library Time: " << (endtime - starttime) << " s\n";

    starttime = omp_get_wtime();
    
    #pragma omp parallel for schedule (static)
    for (int t = 0; t < numberofpoints; t++) 
	{
	    double radians = 2.0 * t * M_PI / numberofpoints;
	    
	    parallelmathlibrary[t] = {static_cast<int>(195 * cos(radians) + screenwidth / 2), static_cast<int>(195 * sin(radians) + screenheight / 2)};
	}
    
    endtime = omp_get_wtime();
    
    cout << "Parallel Math Library Time: " << (endtime - starttime) << " s\n";
    
    vector<SDL_Point> taylorseries(numberofpoints);
    vector<SDL_Point> paralleltaylorseries(numberofpoints);

    starttime = omp_get_wtime();
    
    for (int t = 0; t < numberofpoints; t++) 
    {
        double radians = 2.0 * t * M_PI / numberofpoints;
        
        taylorseries[t] = {static_cast<int>(190 * taylorseriescosine(radians) + screenwidth / 2), static_cast<int>(190 * taylorseriessine(radians) + screenheight / 2)};
    }
    
    endtime = omp_get_wtime();
    
    cout << "Serial Taylor Series Time: " << (endtime - starttime) << " s\n";

    starttime = omp_get_wtime();
    
    #pragma omp parallel for schedule (static)
    for (int t = 0; t < numberofpoints; t++) 
	{
	    double radians = 2.0 * t * M_PI / numberofpoints;
	    
	    paralleltaylorseries[t] = {static_cast<int>(185 * paralleltaylorseriescosine(radians) + screenwidth / 2), static_cast<int>(185 * paralleltaylorseriessine(radians) + screenheight / 2)};
	}
    
    endtime = omp_get_wtime();
    cout << "Parallel Taylor Time: " << (endtime - starttime) << " s\n";
    
	cout << "\nColor Coding:\nWhite - Serial Math Library \nRed - Parallel Math Library \nGreen - Serial Taylor Series \nBlue - Parallel Taylor Series \n\n";
    
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("OpenMP Circle", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, screenwidth, screenheight, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    bool running = true;
    SDL_Event event;
    while (running) 
    {
        while (SDL_PollEvent(&event)) 
        {
            if (event.type == SDL_QUIT) 
            {
                running = false;
            }
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); 
        
        for (const auto& point : mathlibrary) 
        {
        	SDL_RenderDrawPoint(renderer, point.x, point.y);
        }

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); 
        
        for (const auto& point : parallelmathlibrary) 
        {
        	SDL_RenderDrawPoint(renderer, point.x, point.y);
        }

        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255); 
        
        for (const auto& point : taylorseries) 
        {
        	SDL_RenderDrawPoint(renderer, point.x, point.y);
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); 
        
        for (const auto& point : paralleltaylorseries)
        {
        	SDL_RenderDrawPoint(renderer, point.x, point.y);
        }

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

