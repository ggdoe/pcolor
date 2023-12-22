#pragma once
#define SDL_MAIN_HANDLED
#include <stdbool.h>
#include <SDL2/SDL.h>

static bool initialized = false;
SDL_Window * window;
SDL_Renderer * renderer;

static void init_window(int width, int height);
static void close_window();

struct show_state {
    SDL_Event event;
    SDL_Texture *texture;
    SDL_Rect render_rect;
    SDL_Point click_zone_start;
    SDL_Point click_zone_min;
    SDL_Point click_zone_max;
    int pic_w;
    int pic_h;
    double zoom;
    bool redraw;
    bool draw_long_click_zone;
    bool quit;
};

inline static 
void reset_rect(struct show_state *ss)
{
    int w, h;
    SDL_GetWindowSize(window, &w, &h);

    if(w < h * ss->pic_w / ss->pic_h)
    {
        ss->render_rect.x = 0;
        ss->render_rect.y = (h - w * ss->pic_h / ss->pic_w)/2;
        ss->render_rect.w = w;
        ss->render_rect.h = w * ss->pic_h / ss->pic_w;
    }
    else
    {
        ss->render_rect.x = (w - h * ss->pic_w / ss->pic_h)/2;
        ss->render_rect.y = 0;
        ss->render_rect.w = h * ss->pic_w / ss->pic_h;
        ss->render_rect.h = h;
    }
}

inline static 
void center_rect(struct show_state *ss, int cx, int cy)
{
    int w, h;
    SDL_GetWindowSize(window, &w, &h);

    ss->render_rect.x = (double)w / 2.0 - ((double)cx + 0.5) * ss->zoom;
    ss->render_rect.y = (double)h / 2.0 - ((double)cy + 0.5) * ss->zoom;

    if(w < h * ss->pic_w / ss->pic_h)
    {
        ss->render_rect.w = ss->zoom * w;
        ss->render_rect.h = ss->zoom * w * ss->pic_h / ss->pic_w;
    }
    else
    {
        ss->render_rect.w = ss->zoom * h * ss->pic_w / ss->pic_h;
        ss->render_rect.h = ss->zoom * h;
    }
}

void manage_event(struct show_state *ss)
{
    switch (ss->event.type)
        {
            case SDL_WINDOWEVENT:
                if(ss->event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                {
                    ss->zoom = 1.0;
                    reset_rect(ss);
                    ss->redraw = true;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                ss->click_zone_start.x = ss->event.button.x;
                ss->click_zone_start.y = ss->event.button.y;
                ss->click_zone_min = ss->click_zone_start;
                ss->click_zone_max = ss->click_zone_start;
                ss->draw_long_click_zone = true;
                break;
            case SDL_MOUSEMOTION:
                if(ss->event.motion.state == SDL_PRESSED)
                {
                    ss->click_zone_min.x = (ss->click_zone_start.x < ss->event.button.x) ? ss->click_zone_start.x : ss->event.button.x;
                    ss->click_zone_min.y = (ss->click_zone_start.y < ss->event.button.y) ? ss->click_zone_start.y : ss->event.button.y;
                    ss->click_zone_max.x = (ss->click_zone_start.x > ss->event.button.x) ? ss->click_zone_start.x : ss->event.button.x;
                    ss->click_zone_max.y = (ss->click_zone_start.y > ss->event.button.y) ? ss->click_zone_start.y : ss->event.button.y;
                    ss->redraw = true;
                }
                SDL_FlushEvent(SDL_MOUSEMOTION);
                break;
            case SDL_MOUSEBUTTONUP:
                {
                    int w, h;
                    SDL_GetWindowSize(window, &w, &h);

                    double w_ratio = (double)w/(ss->click_zone_max.x - ss->click_zone_min.x);
                    double h_ratio = (double)h/(ss->click_zone_max.y - ss->click_zone_min.y);
                    double new_zoom = (w_ratio < h_ratio) ? w_ratio :  h_ratio;

                    int new_center_x =  ((ss->click_zone_min.x + ss->click_zone_max.x)/2 - ss->render_rect.x) / ss->zoom;
                    int new_center_y =  ((ss->click_zone_min.y + ss->click_zone_max.y)/2 - ss->render_rect.y) / ss->zoom;
                    const double limit_zoom = 20.0;
                    if(new_zoom < limit_zoom)
                        ss->zoom *= new_zoom;

                    center_rect(ss, new_center_x, new_center_y);
                    ss->draw_long_click_zone = false;
                    ss->redraw = true;
                }
                break;
            case SDL_MOUSEWHEEL:
                {
                    int w, h, mouse_x, mouse_y;
                    SDL_GetMouseState(&mouse_x, &mouse_y);
                    SDL_GetWindowSize(window, &w, &h);

                    double old_center_x = ((double)w / 2.0 - (double)ss->render_rect.x);
                    double old_center_y = ((double)h / 2.0 - (double)ss->render_rect.y);
                    double mouse_relative_x = (double)(mouse_x - ss->render_rect.x);
                    double mouse_relative_y = (double)(mouse_y - ss->render_rect.y);
                    const double factor = (ss->event.wheel.y > 0) ? 0.9 : 1.0;
                    
                    int new_center_x = (factor * old_center_x + (1.0 - factor) * mouse_relative_x) / ss->zoom;
                    int new_center_y = (factor * old_center_y + (1.0 - factor) * mouse_relative_y) / ss->zoom;
                    
                    const double zoom_factor = 1.1;
                    ss->zoom *= (ss->event.wheel.y > 0) ? zoom_factor : 1.0/zoom_factor;
                    
                    center_rect(ss, new_center_x, new_center_y);
                    ss->redraw = true;
                }
                break;
            case SDL_KEYDOWN:
                switch(ss->event.key.keysym.sym)
                {
                    case SDLK_r: // reset view
                        reset_rect(ss);
                        ss->zoom = 1.0;
                        ss->redraw = true;
                        break;
                    case SDLK_q:
                        ss->quit = true;
                        break;
                }
                break;
            case SDL_QUIT:
                ss->quit = true;
                break;
        }
}

void draw_state(struct show_state *ss)
{
    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, ss->texture, NULL, &ss->render_rect);
    if(ss->draw_long_click_zone)
    {
        SDL_Rect long_click_zone = {
            ss->click_zone_min.x, 
            ss->click_zone_min.y, 
            ss->click_zone_max.x - ss->click_zone_min.x, 
            ss->click_zone_max.y - ss->click_zone_min.y
        };
        
        SDL_SetRenderDrawColor(renderer, 0xAA, 0xAA, 0xAA, 0xFF);
        SDL_RenderDrawRect(renderer, &long_click_zone);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    }
    SDL_RenderPresent(renderer);
}

struct show_state init_show_state(uint32_t *pixels, int width, int height)
{
    struct show_state ss = {0};
    ss.zoom = 1.0;
    ss.pic_w = width;
    ss.pic_h = height;
    ss.redraw = true;
    
    ss.texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, width, height);
    SDL_UpdateTexture(ss.texture, NULL, pixels, width * sizeof(uint32_t));

    reset_rect(&ss);
    return ss;
}

void show(uint32_t *pixels, int width, int height)
{
    if(!initialized){
        init_window(width, height);
        initialized = true;
    }
    
    struct show_state ss = init_show_state(pixels, width, height);
    
    while(!ss.quit){
        SDL_WaitEvent(&ss.event);
        
        manage_event(&ss);
        if(ss.redraw)
        {
            ss.redraw = false;
            draw_state(&ss);
        }
    } // end main loop

    SDL_DestroyTexture(ss.texture);
}

typedef void(*callback_fn)(void* callback_args);
struct custom_keyevent{
    SDL_Keycode key;
    callback_fn callback;
    void* callback_args;
};

void animate(uint32_t *pixels, int width, int height, double target_fps, callback_fn callback, void* callback_args, struct custom_keyevent *custom_event, int number_custom_keyevent)
{
    if(!initialized){
        init_window(width, height);
        initialized = true;
    }

    const double max_delay = (target_fps == 0) ? 0.0 : 1.0 / target_fps;
    char window_title[256];
    Uint64 delay_rename_title = SDL_GetPerformanceFrequency();
    
    struct show_state ss = init_show_state(pixels, width, height);
    
    while(!ss.quit){
        Uint64 start_ticks = SDL_GetPerformanceCounter();
        while(SDL_PollEvent(&ss.event))
        {
            manage_event(&ss);

            // custom key event
            if(ss.event.type == SDL_KEYDOWN){
                for(int i=0; i < number_custom_keyevent; i++)
                    if(ss.event.key.keysym.sym == custom_event[i].key)
                        custom_event[i].callback(custom_event[i].callback_args);
            }
        }
        callback(callback_args);
        SDL_UpdateTexture(ss.texture, NULL, pixels, width * sizeof(uint32_t));
        draw_state(&ss);

        // active waiting
        while((double)(SDL_GetPerformanceCounter() - start_ticks) / (double)SDL_GetPerformanceFrequency() < max_delay );

        // indicate fps in title twice per second
        if(start_ticks - delay_rename_title > SDL_GetPerformanceFrequency() / 2){
            delay_rename_title = SDL_GetPerformanceCounter();
            double fps = (double)SDL_GetPerformanceFrequency() / (double)(SDL_GetPerformanceCounter() - start_ticks);
            snprintf(window_title, sizeof(window_title), "fps: %.3lf\n", fps);
            SDL_SetWindowTitle(window, window_title);
        }
    } // end main loop

    SDL_DestroyTexture(ss.texture);
}

static void init_window(int width, int height)
{
    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("show",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN);
    renderer = SDL_CreateRenderer(window, -1, 0);
    // renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC | SDL_RENDERER_ACCELERATED);

    atexit(close_window);
}

static void close_window()
{
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

