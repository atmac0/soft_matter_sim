#include "field.h"
#include "data_types.h"

/* PRECONDITION  : None*/
/* POSTCONDITION : A field will object will be initialized in memory, with each element of the field itself set to the empty value (-1)*/
/* */
Field::Field()
{  
  for(int32_t i = 0; i < FIELD_WIDTH; i++)
  {
    for(int32_t j = 0; j < FIELD_HEIGHT; j++)
    {
      field[i][j] = -1;
    }
  }

  frame_num = 1;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
int32_t Field::field_to_png()
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
  
  for(int32_t i = 0; i < FIELD_WIDTH; i++)
  {
    for(int32_t j = 0; j < FIELD_HEIGHT; j++)
    {
      if(field[i][j] == -1)
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

  frame_num++;
  return 0;
}

/* PRECONDITION  : A field has been initialized. */
/* POSTCONDITION : An location in the field will be assigned a particle number. If a particle number was already at the location, 
                   a collision will be recorded before the particle number is assigned. */
/* The point pass in does not have to be normalized to the bounds of the field. If it goes over bounds, it will wrap around and 
   simply come back in through the other side.*/
void Field::place_edge_in_field(coord_t point, field_t particle_num, Collisions* collisions)
{

  point = translate_point_to_field(point);
  
  //If there is already an edge at this point in the field, record the location and particles interacting as a collision.
  if(field[point.x][point.y] != -1 && field[point.x][point.y] != particle_num)
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

/* PRECONDITION  : A field has been initialized. */
/* POSTCONDITION : None */
/* Returns the particle number at a specified point. */
field_t Field::get_particle_at(coord_t point)
{
  point = translate_point_to_field(point);
  return field[point.x][point.y];
}

/* PRECONDITION  : A field has been initialized */
/* POSTCONDITION : A location defined by point in the field will be set to the empty value (-1) */
/* */
void Field::clear_particle_at(coord_t point)
{
  point = translate_point_to_field(point);
  field[point.x][point.y] = -1;
}

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* Translates a point to a value within the indexable bounds of the field. If the point goes over or under the bounds, the point
   will be re-assigned to wrap around and re-enter through the other side of the field. */
coord_t Field::translate_point_to_field(coord_t point)
{
  //Check if edge is inside of the bounds of the field. If not, have it wrap around.
  if(point.x < 0)
  {
    point.x = FIELD_WIDTH + point.x;
  }
  else
  {
    point.x = point.x%FIELD_WIDTH;
  }

  if(point.y < 0)
  {
    point.y = FIELD_HEIGHT + point.y;
  }
  else
  {
    point.y = point.y%FIELD_HEIGHT;
  }

  return point;
}
