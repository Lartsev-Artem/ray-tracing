#include "ray_tracing3d_main.h"

int main(int argc, char* argv[])
{
	std::string name_file_settings;

	if (argc > 1)
		name_file_settings = argv[1];
	else
		name_file_settings = "ray_tracing3d_start_settings_file.txt";

	size_t class_file_vtk;  // задает тип файла --- имена скалярных данных на сетке
	std::string name_file_vtk;
	std::string out_file_grid_vtk;
	std::string out_file_ray_vtk;

	// параметры аккретора и аккреционного диска вне расчётной области
	Type center_sphere[3];
	Type radius_sphere;
	Type internal_radius_disk;
	Type external_radius_disk;

	// параметры картинной плоскости
	Type width_plane;  // безразмерная ширина картинной плоскости
	Type height_plane;  // безразмерная высота картинной плоскости
	int pixels_wide;
	int pixels_high;

	size_t number_of_sort_subdomains;
	size_t number_of_global_steps;  // количесвто конечных кадров 

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, out_file_grid_vtk, out_file_ray_vtk,
		center_sphere, radius_sphere, internal_radius_disk, external_radius_disk,
		width_plane, height_plane, pixels_wide, pixels_high, number_of_sort_subdomains, number_of_global_steps)) {

		std::cout << "Error reading the start settings\n";
		return 1;
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// скалярные данные сетки (unstructured_grid)
	vtkDataArray* density;
	vtkDataArray* absorp_coef;
	vtkDataArray* rad_en_loose_rate;

	clock_t start_clock = clock();
	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	clock_t end_clock = clock();
	std::cout << "\n Reading time of the vtk file: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";

	// один шаг --- одно изображение на картинной плоскости
	for (size_t global_step = 0; global_step < number_of_global_steps; ++global_step) {
		std::cout << "Global step number: " << global_step << '\n';
		
		start_clock = clock();
		OneFullTraceStep(global_step, number_of_global_steps, center_sphere, radius_sphere, internal_radius_disk, external_radius_disk,
			width_plane, height_plane, pixels_wide, pixels_high, number_of_sort_subdomains, unstructured_grid,
			density, absorp_coef, rad_en_loose_rate, out_file_grid_vtk, out_file_ray_vtk);
		clock_t end_clock = clock();
		std::cout << "\n Global step time: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";
	}

	return EXIT_SUCCESS;
}




