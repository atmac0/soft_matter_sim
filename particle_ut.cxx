#include "data_types.h"

void particle_ut(Field* field)
{
  field_t par_num        = 0;
  Particle_type par_type = SQUARE;
  double angular_vel     = 0;
  double orient          = 0;
  double x_vel           = 0;
  double y_vel           = 0;
  coord_t cm_coord;
  cm_coord.x = FIELD_WIDTH/2;
  cm_coord.y = FIELD_HEIGHT/2;
  
  Particle test_par(field, par_num, par_type, angular_vel, orient, x_vel, y_vel, cm_coord);

  translate_to_field_ut(Particle* test_par);
  
}

void translate_to_field_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void find_angle_relative_to_field_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void add_edge_location_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void translate_x_by_1_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void translate_x_by_granular_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void translate_y_by_granular_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void rotate_particle_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void propagate_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void resolve_collisions_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void find_linear_momentum_at_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void find_r_cell_speed_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void is_change_in_linear_in_direction_of_cm(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void is_change_in_angular_in_direction_of_cm_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void determine_direction_of_angular_change_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void find_change_in_momentum_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void draw_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";
}

void delete_edges_ut(Particle* test_par)
{
  std::cout << __func__ << ": FAILURE\n";y
}
