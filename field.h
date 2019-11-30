#pragma once

#include "collisions.h"
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

class Field
{
private:
  uint16_t field[FIELD_WIDTH][FIELD_HEIGHT];
  void initialize_field();  
public:
  Field();
  int32_t field_to_png(uint32_t frame_num);
  void place_edge_in_field(coord_t point, uint32_t particle_num, Collisions* collisions);
};
