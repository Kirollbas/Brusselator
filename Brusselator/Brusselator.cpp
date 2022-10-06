#include <SFML\Graphics.hpp>
#include <iostream>
#include <vector>

// differential equations of Brusselator
void Brusselator(double* x_data , double* fx_data, double* k_data){
    double k3_bx = k_data[2] * (x_data[3] * x_data[0]);
    double k2_x2y = k_data[1] * (x_data[0] * x_data[0] * x_data[1]);

    fx_data[0] = k_data[0] * (x_data[2]) + k2_x2y - k3_bx - k_data[3] * (x_data[0]);
    fx_data[1] = k3_bx - k2_x2y;
};

// common runge cutta method
void common_runge_cutta(double* x_data, double* k_data, double h, void (*f)(double*, double*, double*),std::vector<std::vector<double>> current_runge_kutta_method){
    
    int order = current_runge_kutta_method.size() - 1;

    double** fx_tmp_data = new double*[order];

    for (size_t i = 0; i < order; i++)
        fx_tmp_data[i] = new double[2];
    

    for (size_t column = 0; column < order; column++){

        double* x_tmp_data = new double[2];
        x_tmp_data[0] = x_data[0];
        x_tmp_data[1] = x_data[1];

        for (size_t line = 0; line < column; line++)
        {
            x_tmp_data[0] += h * current_runge_kutta_method[column][line] * fx_tmp_data[column][0];    
            x_tmp_data[1] += h * current_runge_kutta_method[column][line] * fx_tmp_data[column][1];   
        }

        f(x_data,fx_tmp_data[column],k_data);

        delete[] x_tmp_data;
    }
    
    for (size_t i = 0; i < order; i++)
    {
        x_data[0] += h * current_runge_kutta_method[order][i] * fx_tmp_data[i][0];
        x_data[1] += h * current_runge_kutta_method[order][i] * fx_tmp_data[i][1];
    }
    

    for (size_t i = 0; i < order; i++)
        delete[] fx_tmp_data[i];
    delete[] fx_tmp_data;
}

int main(){

    // amount of substance
    double* x_data = new double[4];
    x_data[0] = 3; //X
    x_data[1] = 3; //Y
    x_data[2] = 1; //A
    x_data[3] = 3; //B
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
    double current_time;

    // buthcer's table of current method
    std::vector<std::vector<double>> current_runge_kutta_method = {
        {0},
        {1.0 / 3, 0},
        {-1.0 / 3, 1, 0},
        {1, -1, 1, 0},
        {1.0 / 8, 3.0 / 8, 3.0 / 8, 1.0 / 8}
    };
    

    // for (current_time = 0; current_time < Time; current_time += hop){
    //     std::cout << current_time << ' ' << x_data[0] << ' ' << x_data[1] << '\n';
    //     common_runge_cutta(x_data,k_data,hop,Brusselator,current_runge_kutta_method);
    // }
    
    //Brusselator window

    sf::RenderWindow window(sf::VideoMode(200, 200), "Brusselator");
    window.setFramerateLimit(15);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        common_runge_cutta(x_data,k_data,hop,Brusselator,current_runge_kutta_method);

        double x_div_y = x_data[0] / x_data[1];
        std::cout << x_data[0] << "    " << x_data[1] << '\n';
        window.clear(sf::Color(static_cast<int>(255 * x_div_y),0,static_cast<int>(255 * (1 - x_div_y)),255));
        window.display();
    }


    delete[] x_data;
}