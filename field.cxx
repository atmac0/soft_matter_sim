#include "field.h"
#include "data_types.h"

Field::Field()
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

  char input_name[50];
  char output_name[50];

  //std::string input_name  = "10000x10000.png";
  sprintf(input_name, "%dx%d.png", FIELD_WIDTH, FIELD_HEIGHT);
  sprintf(output_name, "frame_%d.png", frame_num);
  
  sf::Image image;

  if (!image.loadFromFile(input_name))
  {
    return -1;
  }

  sf::Color w(255,255,255);
  sf::Color b(0,0,0);
  
  for(int32_t i=0; i<FIELD_WIDTH; i++)
  {
    for(int32_t j=0; j<FIELD_HEIGHT; j++)
    {
      if(field[i][j] == 0)
      {
	image.setPixel(i,j,w);
      }
      else
      {
	image.setPixel(i,j,b);
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
  if(field[point.x][point.y] != 0 && field[point.x][point.y] != particle_num)
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

uint16_t Field::get_particle_at(coord_t point)
{
  return field[point.x%FIELD_WIDTH][point.y%FIELD_HEIGHT];
}

void Field::clear_particle_at(coord_t point)
{
  field[point.x%FIELD_WIDTH][point.y%FIELD_HEIGHT] = 0;
}
