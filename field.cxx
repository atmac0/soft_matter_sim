#include "field.h"

Field::Field()
{
  initialize_field();
}

void Field::initialize_field()
{  
  for(int32_t i=0; i<FIELD_WIDTH; i++)
  {
    for(int32_t j=0; j<FIELD_HEIGHT; j++)
    {
      field[i][j] = 0;
    }
  } 
}

int32_t Field::field_to_png(uint32_t frame_num)
{

  std::string input_name  = "100x100.png";
  std::string output_name = "frame_1.png";
  sf::Image image;
  //output_name += char(frame_num);

  if (!image.loadFromFile(input_name))
  {
    return -1;
  }

  sf::Color p;
  
  for(int32_t i=0; i<FIELD_WIDTH; i++)
  {
    for(int32_t j=0; j<FIELD_HEIGHT; j++)
    {
      if(field[i][j] == 0)
      {
	p = sf::Color(255,255,255);
	image.setPixel(i,j,p);
      }
      else
      {
	p = sf::Color(0,0,0);
	image.setPixel(i,j,p);
      }
    }
  }

  if (!image.saveToFile(output_name))
  {
    return -1;
  }

  return 0;
}

void Field::place_edge_in_field(coord_t point, uint32_t particle_num, Collisions* collisions)
{
  
  //If there is already an edge at this point in the field, record the location and particles interacting as a collision.
  if(field[point.x][point.y] != 0)
  {
    collision_t collision;
    coord_t collision_coord;
    
    collision_coord.x = point.x;
    collision_coord.y = point.y;
    
    collision.location  = collision_coord;
    collision.particle1 = particle_num;
    collision.particle2 = field[point.x][point.y];
    
    collisions->add_collision(collision);
  }
    
  field[point.x][point.y] = particle_num;
}
