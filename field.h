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
  field_t field[FIELD_WIDTH][FIELD_HEIGHT];
  uint32_t frame_num;
public:
  Field();
  int32_t field_to_png();
  void place_edge_in_field(coord_t point, field_t particle_num, Collisions* collisions);
  field_t get_particle_at(coord_t point);
  void  clear_particle_at(coord_t point);
};
