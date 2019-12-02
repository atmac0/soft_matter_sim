#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "data_types.h"
#include <string>
#include <memory>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "field.h"
#include "particle.h"

int32_t main()
{
 
  Field my_field;
  
  uint32_t par_num1  = 1;
  Particle_type par1 = SQUARE;
  double ang_vel1    = M_PI/100;
  double orient1     = M_PI/4;
  double x_vel1      = 10;
  double y_vel1      = 0;
  coord_t cm_coord1  = {300,500};

  uint32_t par_num2  = 2;
  Particle_type par2 = SQUARE;
  double ang_vel2    = -1;
  double orient2     = 0;
  double x_vel2      = 0;
  double y_vel2      = 0;
  coord_t cm_coord2 = {500,500};  

  Particle my_particle1(&my_field, par_num1, par1, ang_vel1, orient1, x_vel1, y_vel1, cm_coord1);
  Particle my_particle2(&my_field, par_num2, par2, ang_vel2, orient2, x_vel2, y_vel2, cm_coord2); 

  for(uint32_t frame_num=1; frame_num<50; frame_num++)
  {
    my_particle1.propagate();
    my_particle2.propagate();
    my_field.field_to_png(frame_num);
  }
  
  return 0;
}

