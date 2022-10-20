//#include <SFML\Graphics.hpp>
#include <iostream>
#include <cmath>

// differential equations of Brusselator
void Brusselator(double* x_data , double* fx_data, double* k_data){
    double k3_bx = k_data[2] * (x_data[3] * x_data[0]);
    double k2_x2y = k_data[1] * (x_data[0] * x_data[0] * x_data[1]);

    fx_data[0] = k_data[0] * (x_data[2]) + k2_x2y - k3_bx - k_data[3] * (x_data[0]);
    fx_data[1] = k3_bx - k2_x2y;
};

// common runge cutta method
void common_runge_cutta(double* x_data, int vector_dim, int vector_variables ,double* k_data, void(*RS_function)(double*, double*, double*), int steps, double* butchers_table, double hop){
    double* tmp_fx_data = new double[vector_variables * steps];
    double* tmp_x_data = new double[vector_dim];
    

    for (int line = 0; line < steps; line++)
    {   

        for (int i = 0; i < vector_dim; i++)
        {
            tmp_x_data[i] = x_data[i];
        }

        for (int column = 0; column < line; column++)
        {
            for (int i = 0; i < vector_variables; i++)
            {
                tmp_x_data[i] += hop * butchers_table[steps * line + column] * tmp_fx_data[steps * column + i];
            }
            
        }

        RS_function(tmp_x_data, tmp_fx_data + (steps * line), k_data);
        
    }

    for (int i = 0; i < vector_variables; i++)
    {
        for (int step = 0; step < steps; step++)
        {
            x_data[i] += hop * butchers_table[steps * steps + step] * tmp_fx_data[vector_variables * step + i];
        }
    }
    

    delete[] tmp_x_data;
    delete[] tmp_fx_data;
}

int main(){

    // amount of substance
    int width = 1;
    int height = 1;
    double** x_data = new double*[width * height];
    x_data[0] = new double[4];
    x_data[0][0] = 3; //X
    x_data[0][1] = 3; //Y
    x_data[0][2] = 1; //A
    x_data[0][3] = 3; //B
    // 

    // k - reaction velocity
    double* k_data = new double[4];
    k_data[0] = 1; //k1
    k_data[1] = 1; //k2
    k_data[2] = 1; //k3
    k_data[3] = 1; //k4
    //

    //hop
    double hop = 0.01;

    //observation time
    int Time = 50;

    //current time
    double current_time = 0;

    // buthcer's table of euler method 
    int euler_steps = 1;

    double euler_method[2 * 1];
    euler_method[2*0 + 0] = 0; 
    euler_method[2*0 + 1] = 1.0; 

    
    // buthcer's table of RK4 method 
    int RK4_steps = 4;
    double RK4_method[5 * 4];

    RK4_method[4*0 + 0] = 0; 
    RK4_method[4*0 + 1] = 0; 
    RK4_method[4*0 + 2] = 0; 
    RK4_method[4*0 + 3] = 0; 

    RK4_method[4*1 + 0] = (1.0 / 2); 
    RK4_method[4*1 + 1] = 0; 
    RK4_method[4*1 + 2] = 0; 
    RK4_method[4*1 + 3] = 0; 

    RK4_method[4*2 + 0] = 0; 
    RK4_method[4*2 + 1] = (1.0 / 2); 
    RK4_method[4*2 + 2] = 0; 
    RK4_method[4*2 + 3] = 0; 

    RK4_method[4*3 + 0] = 0; 
    RK4_method[4*3 + 1] = 0; 
    RK4_method[4*3 + 2] = 1; 
    RK4_method[4*3 + 3] = 0; 

    RK4_method[4*4 + 0] = (1.0 / 6); 
    RK4_method[4*4 + 1] = (1.0 / 3); 
    RK4_method[4*4 + 2] = (1.0 / 3); 
    RK4_method[4*4 + 3] = (1.0 / 6); 

    


    for (current_time = 0; current_time < Time; current_time += hop){
        std::cout << current_time << ' ' << x_data[0][0] << ' ' << x_data[0][1]<< '\n';
        common_runge_cutta(x_data[0],4,2,k_data,Brusselator,euler_steps,euler_method,hop);
    }
    
    //Brusselator window

    // sf::RenderWindow window(sf::VideoMode(200, 200), "Brusselator");
    // window.setFramerateLimit(30);

    // while (window.isOpen())
    // {
    //     sf::Event event;
    //     while (window.pollEvent(event))
    //     {
    //         if (event.type == sf::Event::Closed)
    //             window.close();
    //     }

    //     common_runge_cutta(x_data,k_data,hop,Brusselator,current_runge_kutta_method);
    //     current_time += hop;
    //     double x_div_y = x_data[0] / x_data[1];
    //     std::cout << current_time << "    "<< x_data[0] << "    " << x_data[1] << '\n';
    //     window.clear(sf::Color(static_cast<int>(255 * x_div_y),0,static_cast<int>(255 * (1 - x_div_y)),255));
    //     window.display();
    // }


    for (size_t i = 0; i < width * height; i++)
    {
        delete[] x_data[i];
    }
    delete[] x_data;
}