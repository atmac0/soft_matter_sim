#include "particle.h"
#include "field.h"

uint16_t Particle::get_particle_mass()
{
  if(particle_type == SQUARE)
  {
    return 1;
  }
}

Particle::Particle(Field* field_, uint16_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord)
{
  particle_num       = par_num;
  particle_type      = par_type;

  angular_velocity   = angular_vel;
  orientation        = orient;
  
  x_velocity         = x_vel;
  y_velocity         = y_vel;

  center_mass_coord  = cm_coord;
  field              = field_;

  if(par_type == SQUARE)
  {
    edge_locations.resize(NUM_SIDES_SQUARE*SQUARE_SIDE_LENGTH/CELL_SIZE); //resize edge location vector to hold all possible edges 
  }
}  

/* The un-rotated orientation of the square is assumed to be as if the top and bottom of the square are parallel with the top and bottom
of the field. A line parallel to these references, pointing to the right, is used as the reference for the x-axis.

........................
.....|-----------|......
.....|...........|......
.....|.....---------->x.
.....|...........|......
.....|-----------|......
........................

The drawing begins at the bottom right corner. This location is found by converting from polar coordinates (using the angle of the orientation),
to cartesian coordinates (using the field as a reference frame), where:

x = rcos(theta)
y = rsin(theta)

where r is calulated by incremented values of an x and y value that are stored relative to the square, using the pythagorean theorem.

*/
void Particle::add_edge_location(coord_t location)
{

  
  
  edge_locations.push_back(location);
}

Collisions* Particle::draw_square()
{
  int32_t side_cell_count = SQUARE_SIDE_LENGTH/CELL_SIZE; /*Gives a normalized size of each cell. Aka, if the side length is 20um, and the cell size is 2um,
							the number of cells needed to draw the full edge will be 20um/2um = 10 cells */

  Collisions* collisions_p = new Collisions(SQUARE);
  coord_t coord_field;
  coord_t coord_square;
  
  for(int32_t i=-side_cell_count/2; i<=side_cell_count/2; i++)
  {
    //right side
    coord_square.x = side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);

    //top side
    coord_square.x = i;
    coord_square.y = -side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
    
    //left side
    coord_square.x = -side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
    
    //bottom side
    coord_square.x = i;
    coord_square.y = side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
  }
  return collisions_p;
}

//Takes a coordinate of an edge point relative to the center of mass of the particle, and translates it into coordinates relative to the field
coord_t Particle::translate_to_field(coord_t edge_point)
{

  /*From a 2d matrix transform:
    x' = xcos(theta) - ysin(theta)
    y' = xsin(theta) + ycos(theta)
    where x and y are the initial coordinates, and x' and y' are the rotated coordinates */
  int32_t x_prime = edge_point.x*cos(-orientation) - edge_point.y*sin(-orientation); 
  int32_t y_prime = edge_point.x*sin(-orientation) + edge_point.y*cos(-orientation);

  //Translate from coordinates relative to the center of mass, to relative to the field
  edge_point.x = center_mass_coord.x + x_prime;
  edge_point.y = center_mass_coord.y + y_prime;	 

  //Check if edge is inside of the bounds of the field. If not, have it wrap around.
  if(edge_point.x < 0)
  {
    edge_point.x = FIELD_WIDTH + edge_point.x;
  }
  else
  {
    edge_point.x = edge_point.x%FIELD_WIDTH;
  }

  if(edge_point.y < 0)
  {
    edge_point.y = FIELD_HEIGHT + edge_point.y;
  }
  else
  {
    edge_point.y = edge_point.y%FIELD_HEIGHT;
  }
  
  
  return edge_point;
}

/* draw edges will identify the particle type, then call the appropriate function to draw the particle */
Collisions* Particle::draw_edges()
{
  if(particle_type == SQUARE)
  {
    return draw_square();
  }
}

void Particle::delete_edges()
{
  for(std::vector<coord_t>::iterator it = edge_locations.begin(); it != edge_locations.end(); it++)
  {
    if(field->get_particle_at(*it) == particle_num)
    {
      field->clear_particle_at(*it);
    }
  }

  edge_locations.clear();  
}

void Particle::propagate()
{

  double time = 0; //time spent propogating

  double x_cell_velocity = abs(x_velocity)/CELL_SIZE; //convert x velocity to cells per micro-second
  double y_cell_velocity = abs(y_velocity)/CELL_SIZE; //convert y velocity to cells per micro-second

  double x_cell_time     = 1/x_cell_velocity; //time spent per cell movement
  double y_cell_time     = 1/y_cell_velocity;

  double x_distance_to_travel = y_cell_time * abs(x_velocity);
  double y_distance_to_travel;

  uint32_t x_cells_to_travel;
  uint32_t y_cells_to_travel;
  
  double x_cell_int_as_double; //used to store the integer portion of the x granular position
  double y_cell_int_as_double; //used to store the integer portion of the y granular position

  double x_granular_decimal; //used to store the portion of the x granular position after the decimal
  double y_granular_decimal; //used to store the portion of the y granular position after the decimal
  
  Collisions* collisions;

  /*propogate the particle for the TIME_INCREMENT amout of time. Curret time is kept track of relative to the change in the y coordinate.
   For every cell moved in the y direction, move x_y_ratio in the x direction. The movement in the y direction will occur first, followed
   by a redraw. For every whole integer in x_y_ratio, the particle will be propagated, with a redraw at each cell, in the x direction that
  number of cells. The rollover will be stored in the granular position variable of the particle.*/
  while(time < TIME_INCREMENT)
  {
    /*Check if the particle will travel atleast 1 cell in the y direction in the time remaining. If it does, proceed by incrementing
      the y direction by one, and the x direction by the number it should travel, determined by the x_y_ratio.*/
    
    if((time + y_cell_time) <= TIME_INCREMENT)
    {

      rotate_particle(y_cell_time);
      
      collisions = translate_y_by_1();
      resolve_collisions(collisions);

      /*being a little tricky here. translate_x_by_granular will return either 0 or 1 if the cell overflows. the modf will take the integer
	portion of the double x_distance_to_travel, which will be assigned to x_cells_to_travel. If the decimal remainder overflows the the cell,
	the number of cells to travel will then have 1 added to it.*/
      
      x_cells_to_travel = 0;
      x_cells_to_travel += translate_x_by_granular(modf(x_distance_to_travel, &x_cell_int_as_double));
      x_cells_to_travel += (uint32_t)x_cell_int_as_double;

      for(uint32_t i=0; i<x_cells_to_travel; i++)
      {
      	collisions = translate_x_by_1();
      	resolve_collisions(collisions);
      }
    }
    /*If incrementing once in the y direction would go over the time increment limit, add the granular amount. If the granular amount goes over one,
      propogate in the y direction. Then propogate x by the correct proportion.*/
    else if((time + y_cell_time) > TIME_INCREMENT)
    {
      double time_remaining = ((double)TIME_INCREMENT) - time;

      rotate_particle(time_remaining);
      
      y_distance_to_travel = abs(y_velocity) * time_remaining;
      x_distance_to_travel = abs(x_velocity) * time_remaining;

      if(translate_y_by_granular(y_distance_to_travel) == 1)
      {
    	collisions = translate_y_by_1();
	resolve_collisions(collisions);
      }
      else
      {
	//This is only called if the y velocity is very small. Without this, the particle will not be drawn if it is stationary.
	delete_edges();
	resolve_collisions(draw_edges());
      }

      x_cells_to_travel = 0;
      x_cells_to_travel += translate_x_by_granular(modf(x_distance_to_travel, &x_cell_int_as_double));
      x_cells_to_travel += (uint32_t)x_cell_int_as_double;

      for(uint32_t i=0; i<x_cells_to_travel; i++)
      {
    	collisions = translate_x_by_1();
    	resolve_collisions(collisions);
      }
    }
    time += y_cell_time;
  }
}

//add the rotational displacement (angular velocity * time span) to the orientation. If the orientation goes over 2*PI, take the remainder.
void Particle::rotate_particle(double time_span)
{
  orientation += angular_velocity * time_span;
  orientation = fmod(orientation, M_PI*2);
}

/*Translate the particle by 1 in the direction appropriate to the current
velocity of the particle in the x direction. This fuction deletes the current edges, propagates
the particle in the appropriate direction by 1, then redraws the edges. This
function returns the collisions that occured from propogation.

If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
the bottom, it will reenter thetop, or if it exits the left, it will reenter the right.*/
Collisions* Particle::translate_x_by_1()
{

  delete_edges();
  
  if(x_velocity > 0)
  {
    center_mass_coord.x++;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }
  else
  {
    center_mass_coord.x--;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }

  return draw_edges();
}

/*Translate the particle by 1 in the direction appropriate to the current
velocity of the particle in the y direction. This fuction deletes the current edges, propagates
the particle in the appropriate direction by 1, then redraws the edges. This
function returns the collisions that occured from propogation.

If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
the bottom, it will reenter thetop, or if it exits the left, it will reenter the right.*/
Collisions* Particle::translate_y_by_1()
{
  
  delete_edges();
  
  if(y_velocity > 0)
  {
    center_mass_coord.y++;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }
  else
  {
    center_mass_coord.y--;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }

  return draw_edges();
}

/*Translates a particle by a granular distance. If the granular postion
  does not exceed the bounds of the cell (is not greater than 1 or less than 0), 
the function will return 0, meaning 0 cells should be moved. If it does exceed those 
bounds, the fuction will return 1, meaning the center of mass postion needs to be moved
by 1.*/
uint32_t Particle::translate_y_by_granular(double granularity)
{
  if(y_velocity > 0)
  {
    y_position_granular += granularity;
  }
  else
  {
    y_position_granular -= granularity;
  }

  if(y_position_granular < 0)
  {
    y_position_granular = 1 + y_position_granular;
    return 1;
  }
  if(y_position_granular > 1)
  {
    y_position_granular = y_position_granular - 1;
    return 1;
  }

  return 0;
}

uint32_t Particle::translate_x_by_granular(double granularity)
{
  if(x_velocity > 0)
  {
    x_position_granular += granularity;
  }
  else
  {
    x_position_granular -= granularity;
  }

  if(x_position_granular < 0)
  {
    x_position_granular = 1 + x_position_granular;
    return 1;
  }
  if(x_position_granular > 1)
  {
    x_position_granular = x_position_granular - 1;
    return 1;
  }

  return 0;
}

void Particle::resolve_collisions(Collisions* collisions)
{
  for(uint32_t i = 0; i < collisions->counts(); i++)
  {
    
    std::cout << "Collision at (" << collisions->get_collision_at(i).location.x << "," << collisions->get_collision_at(i).location.y << ") between particle " << collisions->get_collision_at(i).particle1 << " and particle " << collisions->get_collision_at(i).particle2 << "!\n";
  }
}
