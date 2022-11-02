#include <SFML\Graphics.hpp>
#include <iostream>
#include <cmath>
#include <string>

#define WIDTH 100
#define HEIGHT 100
#define VECTOR_DIM 2

#define CURRENT_X(i,j,n)    2 * (i * (WIDTH + 2) + j) + n
#define LEFT_X(i,j,n)       2 * (i * (WIDTH + 2) + (j - 1)) + n
#define RIGHT_X(i,j,n)      2 * (i * (WIDTH + 2) + (j + 1)) + n
#define UP_X(i,j,n)         2 * ((i - 1) * (WIDTH + 2) + j) + n
#define DOWN_X(i,j,n)       2 * ((i + 1) * (WIDTH + 2) + j) + n
#define CURRENT_FX(i,j,n)   2 * ((i - 1) * WIDTH + (j - 1)) + n


struct Brusselator_data
{
public:
    double* x_data;
    double* fx_data;
    double* k_data;

    Brusselator_data(){
        x_data = new double[(WIDTH + 2) * (HEIGHT + 2) * VECTOR_DIM];
        fx_data = new double[WIDTH * HEIGHT * VECTOR_DIM];
        k_data = new double[8];
    }

    ~Brusselator_data(){
        delete[] x_data;
        delete[] fx_data;
        delete[] k_data;
    }
};


// differential equations of Brusselator
void Brusselator(double* x_data , double* fx_data, double* k_data){
    for (size_t i = 1; i <= HEIGHT; i++)
    {
        for (size_t j = 1; j <= WIDTH; j++)
        {
            double k3_Bx = k_data[2] * (k_data[5] * x_data[CURRENT_X(i,j,0)]); // B + X -> Y + (Nothing)
            double k2_x2y = k_data[1] * (x_data[CURRENT_X(i,j,0)] * x_data[CURRENT_X(i,j,0)] * x_data[CURRENT_X(i,j,1)]); // 2X + Y -> 3X
            double k1_A = k_data[0] * (k_data[4]); // A -> X
            double k4_x = k_data[3] * (x_data[CURRENT_X(i,j,0)]); // X -> (nothing)
            double x_diffusion = k_data[6] * (x_data[LEFT_X(i,j,0)]/4 + x_data[RIGHT_X(i,j,0)]/4 + x_data[UP_X(i,j,0)]/4 + x_data[DOWN_X(i,j,0)]/4 - x_data[CURRENT_X(i,j,0)]); // Diffusion of X
            double y_diffusion = k_data[7] * (x_data[LEFT_X(i,j,1)]/4 + x_data[RIGHT_X(i,j,1)]/4 + x_data[UP_X(i,j,1)]/4 + x_data[DOWN_X(i,j,1)]/4  - x_data[CURRENT_X(i,j,1)]); // Diffusion of Y

            fx_data[CURRENT_FX(i,j,0)] = k1_A + k2_x2y - k3_Bx - k4_x + x_diffusion; // X_FX
            fx_data[CURRENT_FX(i,j,1)] = k3_Bx - k2_x2y + y_diffusion; // Y_FX
        }
        
    }
};


void euler_method(double* x_data, int vec_dim, double* fx_data, double* k_data, void(*RS_function)(double*, double*, double*), double hop){
    RS_function(x_data,fx_data,k_data);


    for (size_t i = 1; i <= HEIGHT; i++)
    {
        for (size_t j = 1; j <= WIDTH; j++)
        {
            for (size_t n = 0; n < vec_dim; n++)
            {
                x_data[CURRENT_X(i,j,n)] += hop * fx_data[CURRENT_FX(i,j,n)];
            }
            
        }
        
    }

    
};

int main(){

    Brusselator_data data;


    for (int i = 0; i < HEIGHT + 2; i++)
    {
        for (int j = 0; j < WIDTH + 2; j++)
        {
            data.x_data[CURRENT_X(i,j,0)] = (rand() % 5) * (rand() % 2) * (rand() % 2);
            data.x_data[CURRENT_X(i,j,1)] = (rand() % 5) * (rand() % 2) * (rand() % 2); 
        }
        
    }
    

    for (int i = 1; i <= WIDTH; i++){

        for (int j = 1; j <= HEIGHT; j++){
            data.x_data[CURRENT_X(i,j,0)] = rand() % 5; //X
            data.x_data[CURRENT_X(i,j,1)] = rand() % 5; //Y
        }
    }


    // k - reaction velocity
    data.k_data[0] = 1; //k1
    data.k_data[1] = 1; //k2
    data.k_data[2] = 1; //k3
    data.k_data[3] = 1; //k4
    data.k_data[4] = 1; //A
    data.k_data[5] = 3; //B
    data.k_data[6] = 0.3; // velocity of X diffusion 
    data.k_data[7] = 0.03; // velocity of Y diffusion 

    //hop
    double hop = 0.05;

    //current time
    double current_time = 0;
    int iterations = 0;

    // sqare size
    int pixel_size_x = 5;
    int pixel_size_y = 5;

    // squares array
    std::vector<std::vector<sf::RectangleShape>> rectangle_shapes;
    rectangle_shapes.resize(WIDTH);

    for (size_t window_width = 0; window_width < WIDTH; window_width++){
        rectangle_shapes[window_width].resize(HEIGHT);

        for (size_t window_height = 0; window_height < HEIGHT; window_height++){
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
    counters.setPosition(0, pixel_size_y * HEIGHT + 5);
    std::string counters_text;

    //Brusselator window

    sf::RenderWindow window(sf::VideoMode(pixel_size_x * WIDTH, pixel_size_y * HEIGHT + 40), "Brusselator");
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

        euler_method(data.x_data,VECTOR_DIM,data.fx_data,data.k_data,Brusselator,hop);
        
        for (size_t i = 1; i <= HEIGHT; i++)
        {
            for (size_t j = 1; j <= WIDTH; j++)
            {   
                double x_div_y = (data.x_data[CURRENT_X(i,j,0)] / data.x_data[CURRENT_X(i,j,1)]) > 1 ? 1 : data.x_data[CURRENT_X(i,j,0)] / data.x_data[CURRENT_X(i,j,1)];
                rectangle_shapes[j - 1][i - 1].setFillColor(sf::Color(static_cast<int>(255 * x_div_y), 0, static_cast<int>(255*(1-x_div_y))));
                window.draw(rectangle_shapes[j - 1][i - 1]);
            }
            
        }

        current_time += hop;
        iterations += 1;

        counters_text = "Time: " + std::to_string(current_time) + "   Iterations: " + std::to_string(iterations);
        counters.setString(counters_text);

        window.draw(counters);
        
        window.display();
    }

}