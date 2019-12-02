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
  char output_name[50];
  sprintf(output_name, "frame_%d.png", frame_num);
  
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

  int32_t field_x; 
  int32_t field_y; 

  if(point.x < 0)
  {
    field_x = FIELD_WIDTH + point.x;
  }
  else
  {
    field_x = point.x%FIELD_WIDTH;
  }

  if(point.y < 0)
  {
    field_y = FIELD_HEIGHT + point.y;
  }
  else
  {
    field_y = point.y%FIELD_HEIGHT;
  }
  
  //If there is already an edge at this point in the field, record the location and particles interacting as a collision.
  if(field[field_x][field_y] != 0 && field[field_x][field_y] != particle_num)
  {
    collision_t collision;
    coord_t collision_coord;
    
    collision_coord.x = field_x;
    collision_coord.y = field_y;
    
    collision.location  = collision_coord;
    collision.particle1 = particle_num;
    collision.particle2 = field[field_x][field_y];
    
    collisions->add_collision(collision);
  }
    
  field[field_x][field_y] = particle_num;
}

uint16_t Field::get_particle_at(coord_t point)
{
  return field[point.x%FIELD_WIDTH][point.y%FIELD_HEIGHT];
}

void Field::clear_particle_at(coord_t point)
{
  field[point.x%FIELD_WIDTH][point.y%FIELD_HEIGHT] = 0;
}
