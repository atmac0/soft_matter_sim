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
  Particle my_particle1(&my_field, 1, SQUARE, M_PI/4, 0, 0, 0, cm_coord1);
  Particle my_particle2(&my_field, 2, SQUARE, 0, 0, 0, 0, cm_coord2); 
  
  Collisions collisions = my_particle1.draw_edges();
  std::cout << collisions.get_collision_count() << "\n";
  collisions = my_particle2.draw_edges();
  std::cout << collisions.get_collision_count() << "\n";
  uint32_t frame_num = 1;
  my_field.field_to_png(frame_num);
  
  return 0;
}

