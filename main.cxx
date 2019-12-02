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
  
  coord_t cm_coord1 = {50,50};
  coord_t cm_coord2 = {45,45};
  
  Particle my_particle1(&my_field, 1, SQUARE, M_PI/100, M_PI/4, 10, 5, cm_coord1);
  //Particle my_particle2(&my_field, 2, SQUARE, 0, 0, 0, 0, cm_coord2); 

  for(uint32_t frame_num=1; frame_num<50; frame_num++)
  {
    my_particle1.propagate();
    my_field.field_to_png(frame_num);
  }
  
  return 0;
}

