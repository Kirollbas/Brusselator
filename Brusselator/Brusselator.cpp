#include <SFML\Graphics.hpp>
#include <iostream>
#include <cmath>
#include <string>



// differential equations of Brusselator
void Brusselator(double* x_data , double* fx_data, double* k_data, double* x_data_left, double* x_data_right, double* x_data_up, double* x_data_down){
    double k3_Bx = k_data[2] * (k_data[5] * x_data[0]); // B + X -> Y + (Nothing)
    double k2_x2y = k_data[1] * (x_data[0] * x_data[0] * x_data[1]); // 2X + Y -> 3X
    double k1_A = k_data[0] * (k_data[4]); // A -> X
    double k4_x = k_data[3] * (x_data[0]); // X -> (nothing)
    double x_diffusion = k_data[6] * (x_data_down[0]/4 + x_data_left[0]/4 + x_data_right[0]/4 + x_data_up[0]/4 - x_data[0]); // Diffusion of X
    double y_diffusion = k_data[7] * (x_data_down[1]/4 + x_data_left[1]/4 + x_data_right[1]/4 + x_data_up[1]/4 - x_data[1]); // Diffusion of Y

    fx_data[0] = k1_A + k2_x2y - k3_Bx - k4_x + x_diffusion;
    fx_data[1] = k3_Bx - k2_x2y + y_diffusion;
};

// common runge cutta method
void common_runge_cutta(double* x_data, int vector_dim ,double* k_data, void(*RS_function)(double*, double*, double*, double*, double*, double*, double*), int steps, double* butchers_table, double hop, double* x_data_left, double* x_data_right, double* x_data_up, double* x_data_down){
    double* tmp_fx_data = new double[vector_dim * steps];
    double* tmp_x_data = new double[vector_dim];
    

    for (int line = 0; line < steps; line++)
    {   

        for (int i = 0; i < vector_dim; i++)
        {
            tmp_x_data[i] = x_data[i];
        }

        for (int column = 0; column < line; column++)
        {
            for (int i = 0; i < vector_dim; i++)
            {
                tmp_x_data[i] += hop * butchers_table[steps * line + column] * tmp_fx_data[steps * column + i];
            }
            
        }

        RS_function(tmp_x_data, tmp_fx_data + (steps * line), k_data, x_data_left, x_data_right, x_data_up, x_data_down);
        
    }

    for (int i = 0; i < vector_dim; i++)
    {
        for (int step = 0; step < steps; step++)
        {
            x_data[i] += hop * butchers_table[steps * steps + step] * tmp_fx_data[vector_dim * step + i];
        }
    }
    

    delete[] tmp_x_data;
    delete[] tmp_fx_data;
}

int main(){

    // amount of substance
    int width = 100;
    int height = 100;

    double** x_data = new double*[(width + 2) * (height + 2)];

    for (int i = 0; i < width + 2; i++)
    {
        for (int j = 0; j < height + 2; j++)
        {
            x_data[i * (width + 2) + j] = new double[2];
            x_data[i * (width + 2) + j][0] = (rand() % 5) * (rand() % 2) * (rand() % 2);
            x_data[i * (width + 2) + j][1] = (rand() % 5) * (rand() % 2) * (rand() % 2); 
        }
        
    }
    

    for (int w = 1; w <= width; w++){

        for (int h = 1; h <= height; h++){
            x_data[w * (width + 2) + h][0] = rand() % 5; //X
            x_data[w * (width + 2) + h][1] = rand() % 5; //Y
        }
    }


    // k - reaction velocity
    double* k_data = new double[8];
    k_data[0] = 1; //k1
    k_data[1] = 1; //k2
    k_data[2] = 1; //k3
    k_data[3] = 1; //k4
    k_data[4] = 1; //A
    k_data[5] = 3; //B
    k_data[6] = 0.2; // velocity of X diffusion 
    k_data[7] = 0.02; // velocity of Y diffusion 

    //hop
    double hop = 0.05;

    //observation time
    int Time = 50;

    //current time
    double current_time = 0;
    int iterations = 0;

    // buthcer's table of euler method 
    int euler_steps = 1;

    double euler_method[2 * 1];
    euler_method[2*0 + 0] = 0; 
    euler_method[2*0 + 1] = 1.0; 



    // sqare size
    int pixel_size_x = 5;
    int pixel_size_y = 5;

    // squares array
    std::vector<std::vector<sf::RectangleShape>> rectangle_shapes;
    rectangle_shapes.resize(width);

    for (size_t window_width = 0; window_width < width; window_width++){
        rectangle_shapes[window_width].resize(height);

        for (size_t window_height = 0; window_height < height; window_height++){
            rectangle_shapes[window_width][window_height].setSize (sf::Vector2f(pixel_size_x,pixel_size_y));
            rectangle_shapes[window_width][window_height].setPosition(pixel_size_x * window_width, pixel_size_y * window_height);
        }
    }
    
    //text of counters;
    sf::Font font;
    font.loadFromFile("calibri_bold.ttf");
    sf::Text counters;
    counters.setFont(font);
    counters.setFillColor(sf::Color(0,255,255));
    counters.setPosition(0, pixel_size_y * height + 5);
    std::string counters_text;

    //Brusselator window

    sf::RenderWindow window(sf::VideoMode(pixel_size_x * width, pixel_size_y * height + 40), "Brusselator");
    window.setFramerateLimit(45);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();

        for (int window_width = 1; window_width <= width; window_width++)
        {
            for (int window_height = 1; window_height <= height; window_height++)
            {
                common_runge_cutta(x_data[window_width * (width + 2) + window_height],2,k_data,Brusselator,euler_steps,euler_method,hop,x_data[(window_width - 1)  * (width + 2) + window_height], x_data[(window_width + 1)  * (width + 2) + window_height], x_data[window_width  * (width + 2) + window_height - 1], x_data[window_width  * (width + 2) + window_height + 1]);
                double x_div_y = (x_data[window_width * (width + 2) + window_height][0] / x_data[window_width * (width + 2) + window_height][1]) > 1 ? 1 : x_data[window_width * (width + 2) + window_height][0] / x_data[window_width * (width + 2) + window_height][1];

                rectangle_shapes[window_width - 1][window_height - 1].setFillColor(sf::Color(static_cast<int>(255 * x_div_y), 0, static_cast<int>(255*(1-x_div_y))));
                window.draw(rectangle_shapes[window_width - 1][window_height - 1]);

            }
        }
        current_time += hop;
        iterations += 1;

        counters_text = "Time: " + std::to_string(current_time) + "   Iterations: " + std::to_string(iterations);
        counters.setString(counters_text);

        window.draw(counters);
        
        window.display();
    }


    for (size_t i = 0; i < width * height; i++)
    {
        delete[] x_data[i];
    }
    delete[] x_data;
}