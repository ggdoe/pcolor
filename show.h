#define SDL_MAIN_HANDLED
#include <stdbool.h>
#include <SDL2/SDL.h>

static bool initialized = false;
SDL_Window * window;
SDL_Renderer * renderer;

static void init_window(int width, int height);
static void close_window();

inline static 
void center_rect(SDL_Rect *rect, int pic_w, int pic_h, int center_x, int center_y, double zoom)
{
    int w, h;
    SDL_GetWindowSize(window, &w, &h);

    // printf("w: %d\th: %d\tpw: %d\tph: %d\tcx: %d\tcy: %d\tzoom: %lf\n", w, h, pic_w, pic_h, center_x, center_y, zoom);

    w *= zoom;
    h *= zoom;

    if(w < h * pic_w / pic_h)
    {
        rect->x = 0;
        rect->y = (h - w * pic_h / pic_w)/2;
        rect->w = w;
        rect->h = w * pic_h / pic_w;
    }
    else
    {
        rect->x = (w - h * pic_w / pic_h)/2;
        rect->y = 0;
        rect->w = h * pic_w / pic_h;
        rect->h = h;
    }

    // TODO: centrer la cam quand on zoom
    rect->x += center_x;
    rect->y += center_y;

    SDL_RenderClear(renderer);
}

void show(uint32_t *pixels, int width, int height)
{
    if(!initialized){
        init_window(width, height);
        initialized = true;
    }
    
    SDL_Rect render_rect;
    SDL_Texture *texture = SDL_CreateTexture(renderer,
            SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, width, height);

    int center_x = 0;
    int center_y = 0;
    double zoom = 1.0;
    const double zoom_factor = 1.1;

    center_rect(&render_rect, width, height, center_x, center_y, zoom);
    SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
    SDL_RenderCopy(renderer, texture, NULL, &render_rect);
    SDL_RenderPresent(renderer);

    SDL_Event event;
    while(SDL_WaitEvent(&event)){
        bool redraw = false;
        switch (event.type)
        {
            case SDL_QUIT:
                goto quit;
            case SDL_WINDOWEVENT:
                if(event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                    redraw = true;
                break;
            case SDL_MOUSEBUTTONUP:
                {
                    int w, h;
                    SDL_GetWindowSize(window, &w, &h);
                    center_x += w/2 - event.button.x;
                    center_y += h/2 - event.button.y;
                    redraw = true;
                }
                break;
            case SDL_MOUSEWHEEL:
                zoom *= (event.wheel.y > 0) ? zoom_factor : 1.0/zoom_factor;
                redraw = true;
                break;
            case SDL_KEYDOWN:
                switch(event.key.keysym.sym)
                {
                    case SDLK_r:
                        center_x = 0;
                        center_y = 0;
                        zoom = 1.0;
                        redraw = true;
                        break;
                    case SDLK_q:
                        goto quit;
                }
                break;
        }
        if(redraw)
        {
            center_rect(&render_rect, width, height, center_x, center_y, zoom);
            SDL_RenderCopy(renderer, texture, NULL, &render_rect);
            SDL_RenderPresent(renderer);
        }
    }
    quit:

    SDL_DestroyTexture(texture);
}

static void init_window(int width, int height)
{
    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("show",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN);
    renderer = SDL_CreateRenderer(window, -1, 0);

    atexit(close_window);
}

static void close_window()
{
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

