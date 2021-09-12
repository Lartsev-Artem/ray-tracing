#include "ray_tracing3d_main.h"
#include "ray_tracing3d_vtk_resourse.h"


// ������������� 3d ����� �� 2d � ��������� ������� local_basis � �������� start_ray
size_t Build2dPlane(const Type* start_ray, Type**& local_basis, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkSmartPointer<vtkPoints>& points_of_plane,
	vtkSmartPointer<vtkUnstructuredGrid>& plane_grid) {

	vtkSmartPointer<vtkCellArray> triangle_cell_array = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> unstructured_grid_points = unstructured_grid->GetPoints();

	Type old_point[3]; // ����� �� �������� �����
	Type new_point[3]; // ����� �� ��������

	// ���������� ��� ����� �� ���������(2d)
	for (size_t i = 0; i < unstructured_grid_points->GetNumberOfPoints(); ++i) {

		unstructured_grid_points->GetPoint(i, old_point);
		Make2dPoint(start_ray, local_basis, old_point, new_point);
		points_of_plane->InsertNextPoint(new_point);
	}

	int size = unstructured_grid->GetNumberOfCells();

	// ������� ����� �� ���������
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < 4; ++j)  // 4 --- ����� ������ 
			triangle_cell_array->InsertNextCell(unstructured_grid->GetCell(i)->GetFace(j));

	plane_grid->SetPoints(points_of_plane); // ������ ����� �� ���������
	plane_grid->SetCells(VTK_TRIANGLE, triangle_cell_array); // ������ ����� �������

	return 0;
}

size_t CheckIntersectionWithGeometryObjects(Type* normal_to_picture_plane, const Type* start_ray, const Type* center_sphere, const Type radius_sphere,
	const Type internal_radius_disk, const Type external_radius_disk) {

	// ������� ����� ������������ (3 �������, 3 ����������)
	std::vector<Type*> cur_face_of_triangle;
	cur_face_of_triangle.resize(3);
	for (size_t i = 0; i < 3; ++i)
		cur_face_of_triangle[i] = new Type[3];

	auto DeleteCurFace{ [&cur_face_of_triangle]() {
			for (size_t i = 0; i < 3; ++i)
		delete[] cur_face_of_triangle[i];
	} };

	std::vector<Type> res = { 0,0,0 };  // ����� == ��������� ���������� ������� � ��������

	//����������� � ������/������
	{
		// ����� �������� ��������� �����
		cur_face_of_triangle[0][0] = 1;
		cur_face_of_triangle[0][1] = 0;
		cur_face_of_triangle[0][2] = 0;

		cur_face_of_triangle[1][0] = 0;
		cur_face_of_triangle[1][1] = 0.9928768384869221;
		cur_face_of_triangle[1][2] = 0.11914522061843064;

		cur_face_of_triangle[2][0] = 2;
		cur_face_of_triangle[2][1] = 0;
		cur_face_of_triangle[2][2] = 0;   // Wolfram

		// ���������� ���� � ���������� �����
		IntersectionWithPlane(cur_face_of_triangle, start_ray, normal_to_picture_plane, res);

		Type local_intersection[3] = { 0,0,0 }; // ����� ����������� � ��������� ����������� ���������

		// �������� �������, �������� ������ ������������� ����� ��� ��������� ���������
		Type vec_1[3] = { 1, 0, 0 };
		Type vec_2[3] = { 0, -0.992877, -0.119145 };  // Wolfram

		for (size_t k = 0; k < 3; k++) { // � ���������
			local_intersection[0] += (res[k]) * vec_1[k];
			local_intersection[1] += (res[k]) * vec_2[k];
		}

		Type* n = normal_to_picture_plane;  // ��� ���������� ������
		Type A = pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2);

		// ����������� ��������� �� ��������� ����������� ����� � ������
		Type radical = 4 * pow((n[0] * (start_ray[0] - center_sphere[0]) + n[1] * (start_ray[1] - center_sphere[1]) + n[2] * (start_ray[2] - center_sphere[2])), 2) -
			4 * A * (1 - 2 * start_ray[0] + pow(start_ray[0], 2) + pow(start_ray[1], 2) + pow(start_ray[2], 2) - pow(radius_sphere, 2));

		if (radical >= 0) {// ���� ����������� �� ������

			Type t = (n[0] - n[0] * start_ray[0] - n[1] * start_ray[1] - n[2] * start_ray[2] - sqrt(radical) / 2) / A;
			Type in_sphere[3] = { n[0] * t + start_ray[0], n[1] * t + start_ray[1], n[2] * t + start_ray[2] };  // ����� �� �����

			// �� ���������� ��������� ����� (� ��������� ����������� (��� �������))
			if (pow(local_intersection[0] - center_sphere[0], 2) + pow(local_intersection[1] - center_sphere[1], 2) <= pow(internal_radius_disk, 2) ||
				pow(local_intersection[0] - center_sphere[0], 2) + pow(local_intersection[1] - center_sphere[1], 2) >= pow(external_radius_disk, 2)) {
				DeleteCurFace();
				return 2;
			}

			// � ��� ��� ����������� ������?
			Type distance_to_spehere = Distance(start_ray[0], start_ray[1], start_ray[2], in_sphere[0], in_sphere[1], in_sphere[2]);
			Type distance_to_plane = Distance(start_ray[0], start_ray[1], start_ray[2], res[0], res[1], res[2]);

			if (distance_to_spehere > distance_to_plane) {
				DeleteCurFace();
				return 1; // ����
			}
			else {
				DeleteCurFace();
				return 2; // �����
			}
		}
		else {  //��� ����������� �� ������ �� ������ �����
			Type square_of_distance = pow(local_intersection[0] - 1, 2) + pow(local_intersection[1], 2);
			if (square_of_distance < pow(external_radius_disk, 2) && (square_of_distance > pow(internal_radius_disk, 2))) {
				DeleteCurFace();
				return 1;
			}
		}
	}

	DeleteCurFace();
	return 0;  // ��� ����������� 
}

size_t FindNumberOfRoots(std::vector<Type>& roots, const size_t n, const  Type a, const Type b, const Type* start, Type*& direction, Type(*f)(const Type*, Type*&, Type)) {
	/* ����� ����������� ������ � ���� �� */

	Type step = (b - a) / (n - 1);
	size_t size = 0;

	Type fi_1 = f(start, direction, a + 0 * step);
	Type fi_2;
	for (size_t i = 1; i < n; ++i) {

		fi_2 = f(start, direction, a + (i)*step);
		//���� ������� �� ������ ������� �� ��� �� => ���������
		if (fi_1 * fi_2 < 0) {
			size++;
			roots.push_back((2 * a + (2. * i - 1) * step) / 2.);
		}
		fi_1 = fi_2;
	}

	std::sort(roots.begin(), roots.end());
	return size;
}

size_t FindIdAndPoints(const std::vector<vtkIdType>& set_subdomain_indexes, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid, const vtkSmartPointer<vtkPoints>& points_of_plane,
	const Type* ray_position_on_plane, Type* normal_to_picture_plane, const Type* start_ray, 
	std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id, size_t& number_of_intersections) {

	if (set_subdomain_indexes.size() == 0) {
		return 0;
	}

	// ������� ����� ������������ (3 �������, 3 ����������)
	std::vector<Type*> cur_face_of_triangle;
	cur_face_of_triangle.resize(3);
	for (size_t i = 0; i < 3; ++i)
		cur_face_of_triangle[i] = new Type[3];

	std::vector<Type> res = { 0,0,0 };  // ����� == ��������� ���������� ������� � ��������

	vtkIdType cur_index;
	size_t size_subdomain = set_subdomain_indexes.size();
	vtkSmartPointer<vtkIdList> idp =
		vtkSmartPointer<vtkIdList>::New();  // �������� ����������

	for (vtkIdType j = 0; j < size_subdomain; ++j) {
		cur_index = set_subdomain_indexes[j];
		idp = plane_grid->GetCell(cur_index)->GetPointIds();

		// ����� �������� ������������ 2d
		for (size_t h = 0; h < 3; ++h)
			points_of_plane->GetPoint(idp->GetId(h), cur_face_of_triangle[h]);

		if (InTriangle(cur_face_of_triangle, ray_position_on_plane)/*��������� ����������*/) {

			//  ����� �������� ������������ 3d
			for (size_t h = 0; h < 3; ++h) {
				// cur_index / 4 --- ����� �� ��������� � 4 ���� ������ (�������� ����������� � ���� ���������� ���������)
				// cur_index % 4 --- ����������, �� ������ ������ (0,1,2,3) => %
				unstructured_grid->GetCell(cur_index / 4)->GetFace(cur_index % 4)->GetPoints()->GetPoint(h, cur_face_of_triangle[h]);
			}

			// � ������������
			IntersectionWithPlane(cur_face_of_triangle, start_ray, normal_to_picture_plane, res);  // local_basis[2] --- �����������(������� � ����. ���������)

			points_and_id[number_of_intersections++] = std::make_pair(res, cur_index / 4);  // ����� ����������� � ������ �� �������� �����
		}
	}

	for (size_t i = 0; i < 3; ++i)
		delete[] cur_face_of_triangle[i];

	return 0;
}

/*���������� �������������� ����� ����.*/
Type MakeIllum(vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id,
	size_t& index, bool is_intersection_with_rosh) {

	if (points_and_id.size() == 0) return 0;  // ��� ����������� 
	if (!density || !absorp_coef || !rad_en_loose_rate) return 0;

	Type point_1[3];
	Type point_2[3];
	vtkIdType id_cell;

	Type I = 0;  // ���������
	Type alpha;  // ����������� ����������
	Type Q;  // ������������� ���������
	Type s;  // ds 
	Type optical_thick = 0;  // ���������� �������

	// ����� ������������� ������������ !!!

	for (vtkIdType i = index - 1; i > 0; --i) {
		id_cell = points_and_id[i - 1].second;  // �������� � ��������������
		Q = rad_en_loose_rate->GetTuple1(id_cell);

		if (I < pow(10, -20) && Q < pow(10, -20)) continue;  // ���� ����� ���, �� ��������� ������������

		for (size_t j = 0; j < 3; ++j) {
			point_1[j] = points_and_id[i].first[j];
			point_2[j] = points_and_id[i - 1].first[j];
		}

		s = MakeLength(point_1, point_2);
		alpha = density->GetTuple1(id_cell) * absorp_coef->GetTuple1(id_cell);
		optical_thick += s * alpha;

		// ���������� �����
		if (alpha > pow(10, -15))
			I = I * exp(-alpha * s) + Q * (1 - exp(-alpha * s)) / alpha;
		else
			I = I * (1 - alpha * s) + Q * s;
	}

	if (is_intersection_with_rosh && optical_thick < 2. / 3) return -3;  // ��� ����� "����" ������ ��� ����
	return I;
}

size_t OneFullTraceStep(const size_t global_step, const size_t number_of_global_steps,
	const Type* center_sphere, const Type radius_sphere, const Type internal_radius_disk, const Type external_radius_disk,
	const Type width_plane, const Type height_plane, const int pixels_wide, const int pixels_high, const size_t number_of_sort_subdomains,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate,
	std::string  out_file_grid_vtk, std::string  out_file_ray_vtk) {

	// ��������� ��������� ������ �� ���������� �������.  0.879518 --- ����� ����, 0.3 --- ������ ������ ��������� ���������
	Type start_angle = 2 * PI / number_of_global_steps * global_step;
	const Type start_ray[3] = { cos(start_angle) + 0.879518, sin(start_angle) , 0.3 };  // ����������� ������� ������������ ������ ���� �� ���� start_angle ������������ �������� ���������
	const Type end_ray[3] = { 0.879518, 0, 0 };

	// ������� ���������� ������
	Type* normal_to_picture_plane;
	Type* vector_1;
	Type* vector_2;

	Type** local_basis;

	SetMemory(3, normal_to_picture_plane, vector_1, vector_2);
	SetMemory(3, local_basis);

	SetDirect(start_ray, end_ray, normal_to_picture_plane);
	SetBasis(start_ray, normal_to_picture_plane, vector_1, vector_2);

	for (size_t i = 0; i < 3; ++i) {
		local_basis[0][i] = vector_1[i];
		local_basis[1][i] = vector_2[i];
		local_basis[2][i] = normal_to_picture_plane[i];
	}

	// ��������� ���������
	vtkSmartPointer<vtkPoints> points_of_plane =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnstructuredGrid> plane_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// ������������� ������� �� ���������
	clock_t start_clock = clock();
	Build2dPlane(start_ray, local_basis, unstructured_grid, points_of_plane, plane_grid);
	clock_t end_clock = clock();

#ifdef CHECK_FULL_TIME
		std::cout << "\n Building time of the plane_grid: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";
		std::cout << " Number of cells in plane_grid: " << plane_grid->GetNumberOfCells() << '\n';
#endif

	// ���������� ����� �� ���������
	start_clock = clock();
	std::vector<std::vector<vtkIdType>> set_subdomain_indexes;
	set_subdomain_indexes.resize(number_of_sort_subdomains * number_of_sort_subdomains); // number_of_sort_subdomains --- �� ������� �� �����������

	size_t max_subdomain = Sort2dGrid(number_of_sort_subdomains, width_plane, height_plane, plane_grid, points_of_plane, set_subdomain_indexes);
	end_clock = clock();

#ifdef CHECK_FULL_TIME
		std::cout << "\n Sorting time of cells on the plane: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << '\n';
#endif

	// ����������� � ���������� ��������� ���������
	start_clock = clock();
	vtkSmartPointer<vtkUnstructuredGrid> image_plane_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkUnstructuredGrid> ray_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	RayTracingAndBuildImage(points_of_plane, plane_grid, image_plane_grid, ray_grid, unstructured_grid, density, absorp_coef, rad_en_loose_rate, number_of_sort_subdomains,
		max_subdomain, set_subdomain_indexes, start_ray, local_basis, normal_to_picture_plane, width_plane, height_plane, pixels_wide, pixels_high,
		center_sphere, radius_sphere, internal_radius_disk, external_radius_disk);
	end_clock = clock();

#ifdef CHECK_FULL_TIME
		std::cout << "\n Full trace time : " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << '\n';
#endif

#ifdef DRAW_RAY
		out_file_ray_vtk.pop_back();
		out_file_ray_vtk += std::to_string(global_step);
		WriteFileVTK(out_file_ray_vtk + ".vtk", ray_grid);
#endif

	out_file_grid_vtk.pop_back();
	out_file_grid_vtk += std::to_string(global_step);
	WriteFileVTK(out_file_grid_vtk + ".vtk", image_plane_grid);

	ClearMemory(normal_to_picture_plane, vector_1, vector_2);
	ClearMemory(3, local_basis);

	return 0;
}

size_t RayTracingAndBuildImage(const vtkSmartPointer<vtkPoints>& points_of_plane, const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid,
	vtkSmartPointer<vtkUnstructuredGrid>& image_plane_grid, vtkSmartPointer<vtkUnstructuredGrid>& ray_grid,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate,
	const size_t number_of_sort_subdomains, const size_t max_subdomain, const  std::vector<std::vector<vtkIdType>>& set_subdomain_indexes,
	const Type* start_ray, Type**& local_basis, Type*& normal_to_picture_plane, const Type width_plane, const Type height_plane, const int pixels_wide, const int pixels_high,
	const Type* center_sphere, const Type radius_sphere, const Type internal_radius_disk, const Type external_radius_disk) {

	// ���������� �����
#ifdef DRAW_RAY  
	vtkSmartPointer<vtkPoints> points_line =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> line_array =
		vtkSmartPointer<vtkCellArray>::New();
#endif

	// ���������� ��������� ���������
	vtkSmartPointer<vtkPoints> points_quad =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> quad_array =
		vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkDoubleArray> Illum_array =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkQuad> quad =
		vtkSmartPointer<vtkQuad>::New();


	std::vector<std::pair<std::vector<Type>, vtkIdType>> points_and_id;  // ����: ����� ��������� ���� � ������(3d) ������ � id 
	points_and_id.resize(max_subdomain);


	// ������� � Eigen �������� ������ ��� ������������� ������� ����������.
	Eigen::Matrix3d basis;
	Eigen::Vector3d origin_of_coord = { start_ray[0], start_ray[1], start_ray[2] };
	for (size_t i = 0; i < 3; ++i) {
		basis(i, 0) = local_basis[0][i];  //vector_1[i];
		basis(i, 1) = local_basis[1][i];  //vector_2[i];
		basis(i, 2) = local_basis[2][i];  //normal_to_picture_plane[i];
	}

	const Eigen::Vector3d angle_of_picture_plane = { -(width_plane / 2), -(height_plane / 2), 0 };  // ���� ���������. �� ���� ���������� ���������� ���� ���������	
	Eigen::Vector3d Eigen_cur_ray_start;  // ������ ���� �� ������ ���� (����� �������)

	// �������� ������� (2�� �.�. ����� �� ��������� � ��������� �����������)
	const Eigen::Vector3d e1 = { 1,0,0 };
	const Eigen::Vector3d e2 = { 0,1,0 };

	Type step = width_plane / pixels_wide;  // ������ ������� (���==height_plane/pixels_high;)
	Type cur_ray_start[3];
	Type cur_ray_end[3];
	Type cur_ray_position_on_plane[3];  // �������� �������� ������� --- �������� �������� ���� �� �������� ���������

	size_t quad_number = 0;
	size_t ray_number = 0;

	size_t number_of_intersections = 0;  // ����� ���������� ���� � �����

	std::vector<Type> roots_rosh_equation;
	bool is_intersection_with_rosh = false;

	for (size_t i = 0; i < pixels_wide; ++i)
		for (size_t j = 0; j < pixels_high; ++j)
		{
			roots_rosh_equation.resize(0);
			is_intersection_with_rosh = false;
			number_of_intersections = 0;

			// ������� � ������ ������ �������
			Eigen_cur_ray_start[0] = angle_of_picture_plane(0) + i * step;
			Eigen_cur_ray_start[1] = angle_of_picture_plane(1) + j * step;
			Eigen_cur_ray_start[2] = 0;

			for (size_t h = 0; h < 3; ++h)
				cur_ray_position_on_plane[h] = Eigen_cur_ray_start(h);  // �� ���������

			Eigen_cur_ray_start = basis * Eigen_cur_ray_start + origin_of_coord;  // ������� � 3d

			// ������� ����� � ������ ���� � ������������
			for (size_t h = 0; h < 3; ++h) {
				cur_ray_start[h] = Eigen_cur_ray_start(h);
				cur_ray_end[h] = Eigen_cur_ray_start(h) + 4 * normal_to_picture_plane[h];  // 4 --- ������� ����� ����, �����������, ����� ���� ����� �� �������
			}

			// ���������� ������ ����
#ifdef DRAW_RAY  
			points_line->InsertNextPoint(cur_ray_start);
			points_line->InsertNextPoint(cur_ray_end);
			vtkSmartPointer<vtkLine> line =
				vtkSmartPointer<vtkLine>::New();
			line->GetPointIds()->SetId(0, ray_number++);
			line->GetPointIds()->SetId(1, ray_number++);
			line_array->InsertNextCell(line);
#endif  // DRAW_RAY

			// ���������� ������ �������
			{
				points_quad->InsertNextPoint(cur_ray_position_on_plane[0] - step / 2, cur_ray_position_on_plane[1] - step / 2, 0);
				points_quad->InsertNextPoint(cur_ray_position_on_plane[0] + step / 2, cur_ray_position_on_plane[1] - step / 2, 0);
				points_quad->InsertNextPoint(cur_ray_position_on_plane[0] + step / 2, cur_ray_position_on_plane[1] + step / 2, 0);
				points_quad->InsertNextPoint(cur_ray_position_on_plane[0] - step / 2, cur_ray_position_on_plane[1] + step / 2, 0);


				quad->GetPointIds()->SetId(0, quad_number++);
				quad->GetPointIds()->SetId(1, quad_number++);
				quad->GetPointIds()->SetId(2, quad_number++);
				quad->GetPointIds()->SetId(3, quad_number++); // (��������, ������)
				quad_array->InsertNextCell(quad);
			}

			size_t number_of_roots = FindNumberOfRoots(roots_rosh_equation, 1001, -3, 3, cur_ray_start, normal_to_picture_plane, Rosh); // 1001 --- ����� �����, [-3,3] --- �������� ������
			if (number_of_roots == 6 || number_of_roots == 4 /*4,6 --- ����� � ������������ �����������*/) {

				Type intersection_point[3] = { cur_ray_start[0] + roots_rosh_equation[1] * normal_to_picture_plane[0],
						                   	   cur_ray_start[1] + roots_rosh_equation[1] * normal_to_picture_plane[1],
							                   cur_ray_start[2] + roots_rosh_equation[1] * normal_to_picture_plane[2] };

				if (intersection_point[0] < 0.3)  // 0.3 --- ����� � ���������� ����� ������ � ����������� (����������� � ������������ ������, � �� ���������)
					is_intersection_with_rosh = true;
			}

			// test --- �������� �����������, ��� ������ ����. (�����, ��������, ���� ��� �������, ����� �������)
			int test = CheckIntersectionWithGeometryObjects(normal_to_picture_plane, cur_ray_start, center_sphere, radius_sphere,
				internal_radius_disk, external_radius_disk);
					
			if (!test) {
				// ����� ������ ���������� (� ���� ���������� ��� �������� "�� �������")

				int num_x = (int)((cur_ray_position_on_plane[0] + width_plane / 2) / (width_plane / number_of_sort_subdomains));
				int num_y = (int)((cur_ray_position_on_plane[1] + height_plane / 2) / (height_plane / number_of_sort_subdomains));
				
				test = FindIdAndPoints(set_subdomain_indexes[num_x * number_of_sort_subdomains + num_y], unstructured_grid, plane_grid, points_of_plane,
					cur_ray_position_on_plane, normal_to_picture_plane, cur_ray_start, points_and_id, number_of_intersections);

				SortIntersectionPoints(cur_ray_start, number_of_intersections, points_and_id);
			}

			if (test > 0) {
				// ����������� � �������������� ��������
				Illum_array->InsertNextTuple1(-test);  // - ��� ��������� ParaView
			}
			else {
				Type I = MakeIllum(density, absorp_coef, rad_en_loose_rate, points_and_id, number_of_intersections, is_intersection_with_rosh);
				Illum_array->InsertNextTuple1(I);
			}

		}  // for i,j

	image_plane_grid->SetPoints(points_quad);  // ������ �����
	image_plane_grid->SetCells(VTK_QUAD, quad_array);  // ������ �����
	image_plane_grid->GetCellData()->SetScalars(Illum_array);

#ifdef DRAW_RAY
	ray_grid->SetPoints(points_line);  // ������ �����
	ray_grid->SetCells(VTK_QUAD, line_array);  // ������ �����
#endif

	return 0;
}

size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	/*������ ��������� ����� � ������ � vtkUnstructuredGrid*/

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	switch (class_file_vtk) {
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	}

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}

size_t Sort2dGrid(const size_t N/*����� ����������� �� ������� �����������*/, const Type width_plane, const Type height_plane, const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid,
	const vtkSmartPointer<vtkPoints>& points_of_plane, std::vector<std::vector<vtkIdType>>& set_subdomain_indexes) {
	/*���� ����������: ������� �������� �� ��������������� � ����� ��� ������ ������������� ��������������*/

	size_t max_subdomain = 0;
	int number_of_unordered_cells = plane_grid->GetNumberOfCells();  // ���������� ����� ����� 
	std::vector<vtkIdType> ind(number_of_unordered_cells);  // ��� ������� �����
	for (vtkIdType i = 0; i < number_of_unordered_cells; ++i)
		ind[i] = i;

	Type left_border_of_domain[3] = { -width_plane / 2, -height_plane / 2, 0 };  // ��������� "�����" ������� ���������� (������������� --- ���� ��������� ���������)
	Type w_step = width_plane / N;  // ������ ����������
	Type h_step = height_plane / N;  // ������ ����������

	// ����� ��������� ������ �����. true --> ������ � ��������  ������ �� ��� ��������� ����������
	bool is_nonX = false;
	bool is_nonY = false;

	auto cmpXX{ [points_of_plane, w_step, &is_nonX](const Type* left_border_of_domain, const vtkCell* left) {
		/*��������� �� ������ ����������*/
		Type a[3], b[3], c[3];
		vtkIdList* ids = left->PointIds;
		points_of_plane->GetPoint(ids->GetId(0), a);
		points_of_plane->GetPoint(ids->GetId(1), b);
		points_of_plane->GetPoint(ids->GetId(2), c);

		Type A = a[0];
		Type B = b[0];
		Type C = c[0];
		Type X = left_border_of_domain[0];  // "�����" ������� �����������

		// ������ ������� �� ����������� �� ����� ����������
		if (A < X && B < X && C < X) {
			is_nonX = true;
			return false;
		}

		// ������� ����� ������ ���������� (true\false)
		bool eq1 = (X < A&& A < X + w_step);
		bool eq2 = (X < B&& B < X + w_step);
		bool eq3 = (X < C&& C < X + w_step);

		if (eq1 && eq2 && eq3) {
			is_nonX = true;
			return true;
		}

		if (eq1 || eq2 || eq3) {
			return true;
		}

		if (A < X && X < C)
			return true;
		else if (A < X && X < B)
			return true;
		else if (B < X && X < C)
			return true;
		else if (B < X && X < A)
			return true;
		else if (C < X && X < B)
			return true;
		else if (C < X && X < A)
			return true;
		else return false;
} };

	auto cmpYY{ [points_of_plane, h_step, &is_nonY](const Type* left_border_of_domain, const vtkCell* left) {
		/*������ cmpXX, ������ �� ������� �����������*/
		Type a[3], b[3], c[3];
		vtkIdList* ids = left->PointIds;
		points_of_plane->GetPoint(ids->GetId(0), a);
		points_of_plane->GetPoint(ids->GetId(1), b);
		points_of_plane->GetPoint(ids->GetId(2), c);

		Type A = a[1];
		Type B = b[1];
		Type C = c[1];
		Type X = left_border_of_domain[1];

		if (A < X && B < X && C < X) {
			is_nonY = true;
			return false;
		}

		bool eq1 = (X < A&& A < X + h_step);
		bool eq2 = (X < B&& B < X + h_step);
		bool eq3 = (X < C&& C < X + h_step);

		if (eq1 && eq2 && eq3) {
			is_nonY = true;
			return true;
		}

		if (eq1 || eq2 || eq3) {
			return true;
		}

		if (A < X && X < C)
			return true;
		else if (A < X && X < B)
			return true;
		else if (B < X && X < C)
			return true;
		else if (B < X && X < A)
			return true;
		else if (C < X && X < B)
			return true;
		else if (C < X && X < A)
			return true;
		else return false;
} };


	// ��������� �� ����������	
	for (size_t k = 0; k < N; ++k) {
		for (size_t h = 0; h < N; ++h) {
			// ���������� ������������ ����������
			for (size_t i = 0; i < number_of_unordered_cells; ++i) {
				if (cmpXX(left_border_of_domain, plane_grid->GetCell(ind[i])) && cmpYY(left_border_of_domain, plane_grid->GetCell(ind[i])))
					set_subdomain_indexes[k * N + h].push_back(ind[i]);

				// ���� ������ �� ����������� �� ����� ���������� => ��������� �� �� �������
				if (is_nonX && is_nonY) {
					is_nonX = false;
					is_nonY = false;
					std::swap(ind[i--], ind.back());
					--number_of_unordered_cells;
					ind.pop_back();
				}
			}

			left_border_of_domain[1] += h_step;

			if (set_subdomain_indexes[k * N + h].size() > max_subdomain)
				max_subdomain = set_subdomain_indexes[k * N + h].size();
		}
		left_border_of_domain[0] += w_step;
		left_border_of_domain[1] -= N * h_step;
	}
	return max_subdomain;
}

size_t SortIntersectionPoints(const Type* start_ray, size_t& number_of_intersections, std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id) {

	if (number_of_intersections == 0) {
		return 0;
	}

	auto CmpDistancePoint{ [start_ray](const std::pair<std::vector<Type>, vtkIdType >& left, const std::pair<std::vector<Type>, vtkIdType >& right) {
		Type L1 = 0;
		Type L2 = 0;
		for (size_t i = 0; i < 3; ++i) {
			L1 += pow(start_ray[i] - left.first[i], 2);
			L2 += pow(start_ray[i] - right.first[i], 2);
		}
		return L1 < L2;  // �������������� ����������
} };

	std::sort(points_and_id.begin(), points_and_id.begin() + number_of_intersections, CmpDistancePoint);

	// ��������� ������������� �������� ������(1-1, 2-2, 3-3  � �.�) (�.�. �� ����� ���� ����� � 2 ���� ������ (��� ����� �� i ����� � j))
	for (size_t i = 0; i < number_of_intersections - 1; i += 2) {
		if (points_and_id[i].second != points_and_id[i + 1].second)
			swap(points_and_id[i + 1], points_and_id[i + 2]);
	}

	int buf = number_of_intersections - 1;
	number_of_intersections = number_of_intersections / 2 + 1;

	// ��������� ���������� ��������� � ������ (��������� ������� ��������, �� ��� �� ��������� number_of_intersections)
	for (size_t i = 1; i < number_of_intersections - 1; ++i) {
		points_and_id[i] = points_and_id[2 * i];
	}

	points_and_id[number_of_intersections - 1] = points_and_id[buf];
	return 0;
}

size_t WriteFileVTK(const std::string out_file_grid_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	vtkSmartPointer<vtkGenericDataObjectWriter> writer_vtk =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer_vtk->SetFileName(out_file_grid_vtk.c_str());
	writer_vtk->SetInputData(unstructured_grid);
	if (writer_vtk->Write())
		return 0;

	std::cout << "Error writing to the vtk file\n";
	return 1;
}










