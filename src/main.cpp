#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <utility>
#include <tuple>
#include <memory>
#include <vector>
#include <cmath>
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;
using std::pair;
using std::pow;
using std::sqrt;
using std::abs;
#include "io.h"
#include "matrix.h"

#include "MyObject.h"
typedef long double ldouble;
int di[] = { 1, -1, 0, 0 };
int dj[] = { 0, 0, 1, -1 };

void bur(const Image & in, Matrix < uint > & bn) {
	uint n_r = in.n_rows, n_c = in.n_cols;
	uint r, g, b, h;
	h = 175;
	for (uint i = 0; i < n_r; i++)
		for (uint j = 0; j < n_c; j++) {
			tie(r, g, b) = in(i, j);
			if (r > h || g > h || b > h)
				bn(i, j) = 1;
			else
				bn(i, j) = 0;
		}
}

//
// размечивает шестеренку числом - с
void fll(const Matrix< uint > &b, Matrix< uint > &l, int i, int j, int c) {
	int n_r = b.n_rows, n_c = b.n_cols;
	if (i < 0 || i == n_r) return;
	if (j < 0 || j == n_c) return;
	if (l(i, j) || !b(i, j)) return;
	l(i, j) = c;
	fll(b, l, i + 1, j, c);
	fll(b, l, i, j + 1, c);
	fll(b, l, i - 1, j, c);
	fll(b, l, i, j - 1, c);
}

// размечивает все шестеренки 1 2 ... 
int lb(Matrix< uint > & b, Matrix< uint > & l) {
	uint n_r = b.n_rows;
	uint n_c = b.n_cols;
	char c = 0;
	for (uint i = 0; i < n_r; i++)
		for (uint j = 0; j < n_c; j++)
			if (b(i, j) && !l(i, j))
				fll(b, l, i, j, ++c);
	return c;
}
/*
void dfs(uint i, uint j, Matrix< uint > & visit, const Matrix< uint > lbl_Img, uint c) {
	visit(i, j) += 1;
	for (uint k = 0; k < 4; k++) {
		uint i1 = i + di[k];
		uint j1 = j + dj[k];
		if ((i1 < lbl_Img.n_rows && i1 >= 0) && (j1 < lbl_Img.n_cols && j1 >= 0))
			if (lbl_Img(i1, j1) == c && !visit(i1, j1))
				dfs(i1, j1, visit, lbl_Img, c);
	}
}

void define_dot(Matrix< uint > & lbl_Img) {
	int n_r = lbl_Img.n_rows, n_c = lbl_Img.n_cols;
	Matrix< uint > visit(n_r, n_c);
	for (uint i = 0; i < n_r; i++)
		for (uint j = 0; j < n_c; j++)
			visit(i, j) = 0;
	
	for (uint i = 0; i < n_r; i++)
		for (uint j = 0; j < n_c; j++)
			if (!visit(i, j) && !lbl_Img(i, j))
				dfs(i, j, visit, lbl_Img, lbl_Img(i, j));
			
}*/

void print_area(const vector< uint > & v)
{
	for(size_t i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}	
void get_area(vector< uint > & v, uint n, Matrix< uint > img)
{
	v.resize(n+1);
	for(size_t i = 0; i < img.n_rows; i++)
		for(size_t j = 0; j < img.n_cols; j++)
		{
			if (img(i, j))
				v[img(i, j)] += 1;			
		}
	print_area(v);
}

void print_pos(const vector<pair< uint, uint>> & pos)
{
	cout << endl;
	for(size_t i = 0; i < pos.size(); i++)
	{
		cout << i << " (x, y) : (" << pos[i].first << ", ";
		cout << pos[i].second << ")" << endl;
	}
}
void find_centers(vector< uint > & area, Matrix< uint> & img, vector<pair< uint, uint >> & pos)
{
	for(uint i = 0; i < img.n_rows; i++)
		for(uint j = 0; j < img.n_cols; j++)
		{
			if (img(i, j)){
				pos[img(i,j)].first += i;
				pos[img(i,j)].second += j;
			}
		}
	for(size_t i = 1; i < pos.size(); i++)
	{
		pos[i].first = pos[i].first / area[i];
		pos[i].second = pos[i].second / area[i];
	}
	
}/*
void get_radius(vector< int > & area, vector<pair< int, int >> & pos, const Matrix< uint > & img, double & R, double & r)
{
	int num_min = 1;
	for(size_t i = 2; i < area.size(); i++)
		if (num_min > area[i])
			num_min = i;
	pair< int, int > num_gear1;
	//pair< int, int > num_gear2;
	double distance = 1e+9;
	for(size_t i = 1; i < pos.size(); i++){
		if (num_min != i){
			double tmp = sqrt(pow((pos[i].first - pos[num_min].first), 2) -
						pow((pos[i].second - pos[num_min].second), 2));
			if (distance > tmp){
				distance = tmp;
				num_gear1.first = pos[i].first;
				num_gear1.second = pos[i].second;
			}
		}

	for(uint i = 0; i < img.n_rows; i+)
		for(uint j = 0; j < img.n_cols; j+){
			if (img(i, j) == ){
				
		}
}*/
void get_distance(vector< uint > & area, vector<pair< uint, uint >> & pos, double & distance, uint & num_min, uint & num_gear1_img)
{
	uint s;
	if (area.size())
		s = area[1];
	num_min = 1;
	for(size_t i = 2; i < area.size(); i++)
		if (s > area[i]){
			s = area[i];
			num_min = i;
		}
	cout << "\ndot area : " << area[num_min] << endl;
	//pair< int, int > num_gear2;
	
	//distance = sqrt(pow((pos[1].first - pos[num_min].first), 2) -
	//					pow((pos[1].second - pos[num_min].second), 2));
	distance = 1e+9;
	for(size_t i = 1; i < pos.size(); i++){
		if (num_min != i){
			double tmp = sqrt((pos[i].first - pos[num_min].first) * (pos[i].first - pos[num_min].first)+
						(pos[i].second - pos[num_min].second)*(pos[i].second - pos[num_min].second));
			if (distance > tmp){
				distance = tmp;
				num_gear1_img = i;
			}
		}
	}
}

void get_radius(uint x, uint y, Matrix< uint > & lblImg, double & R, double & r, uint num)
{
	double max, min;
	max = 0.0;
	min = sqrt(pow(x, 2) + pow(y, 2));
	for(uint i = 0; i < lblImg.n_rows; i++)
		for(uint j = 0; j < lblImg.n_cols; j++){
			double tmp = sqrt((i-x)*(i-x) + (j-y)*(j-y));
			if (lblImg(i, j) == num && x != i && y != j && tmp > max)
				max = tmp;
			
			if (!lblImg(i, j) && x != i && y != j && tmp < min)
				min = tmp;
		}
	R = max;
	r = min;
	/*cout << endl << endl;
	for(uint i = 0; i < lblImgGear.n_rows; i++){
		for(uint j = 0; j < lblImgGear.n_cols; j++)
			cout << lblImgGear(i, j);
		cout << endl;
	}*/
}

void paste(double R[], double distance, Matrix< uint > & lbl_Img, Matrix< uint > & lbl_ImgGear, uint x, uint y)
{
	if (R[1] + R[2] < distance && R[0] + R[2] > distance)
		cout << "suitable" << endl;
	else{
		cout << "bad gear" << endl;
		return;
	}
	cout << "paste : ";
	cout << x << " " << y << endl;
	uint n_r = lbl_Img.n_rows, n_c = lbl_Img.n_cols; 
	for(uint i = 0; i < lbl_ImgGear.n_rows; i++)
		for(uint j = 0; j < lbl_ImgGear.n_cols; j++){
			if (i + x < n_r && j + y < n_c)
				if (lbl_ImgGear(i, j)){
					lbl_Img(i+x, j+y) = 100;
				}
		}
	
}

tuple <int, vector< shared_ptr < IObject >>, Image, Image >
repaire(Image & in, Image & gear) { // разметка шестеренок|  num - имя.bmp
	Matrix< uint > binImg(in.n_rows, in.n_cols), lbl_Img(in.n_rows, in.n_cols);
	Matrix< uint > binImgGear(gear.n_rows, gear.n_cols), lbl_ImgGear(gear.n_rows, gear.n_cols);
	auto object_array = vector< shared_ptr < IObject >>(); // автоматичекская определение и инициализация
	uint result_idx = 0;
	uint n_r = in.n_rows, n_c = in.n_cols;
	uint nobj = 0, nobjg = 0;
	bur(in, binImg); //бинаризация
	bur(gear, binImgGear);
	for(uint i = 0; i < n_r; i++)
		for(uint j = 0; j < n_c; j++)
			lbl_Img(i, j) = 0;
	for(uint i = 0; i < gear.n_rows; i++)
		for(uint j = 0; j < gear.n_cols; j++)
			lbl_ImgGear(i, j) = 0;
	nobj = lb(binImg, lbl_Img); // разметка
	nobjg = lb(binImgGear, lbl_ImgGear);
	cout << "shesterenok : " << nobj - 1 << endl; // не считая точк
	cout << "gear num : " << nobjg << endl;
	if (nobjg != 1)
		throw "must have only ONE gear!!!";
	vector< uint > area;
	vector< uint > areag;
	get_area(area, nobj, lbl_Img);
	get_area(areag, nobjg, lbl_ImgGear);
	vector<pair <uint, uint>> pos(area.size());
	vector<pair <uint, uint>> posg(areag.size());
	find_centers(area, lbl_Img, pos);
	find_centers(areag, lbl_ImgGear, posg);
	print_pos(pos);
	print_pos(posg);
	cout << "gear center : " << posg[1].first << " " << posg[1].second << endl;
	double distance;
	uint num_dot;// їїїїї їїїїїї їїїїї
	uint num_gear1_img;
	get_distance(area, pos, distance, num_dot, num_gear1_img);
	cout << "\ndot pos : " << pos[num_dot].first << " " << pos[num_dot].second << endl;
	cout << "\ndistance : " << distance << endl;	
	double R_img_gear, r_img_gear;
	double R_img, r_img;
	get_radius(pos[num_gear1_img].first, pos[num_gear1_img].second, lbl_Img, R_img, r_img, num_gear1_img);
	get_radius(posg[1].first, posg[1].second, lbl_ImgGear, R_img_gear, r_img_gear, 1);	
	cout << "\npos gear (img) : " << pos[num_gear1_img].first;
	cout << " " << pos[num_gear1_img].second << endl;
	cout << "R (img) : " <<  R_img << " r (img) : " << r_img << endl;
	cout << "R : " << R_img_gear << " r : " << r_img_gear << endl;
	double radius[4] = {
		R_img,
		r_img,
		R_img_gear,
		r_img_gear
	};
	uint x = abs(pos[num_dot].first - posg[1].first);
	uint y = abs(pos[num_dot].second - posg[1].second);
	paste(radius, distance, lbl_Img, lbl_ImgGear, x, y);
	Image res(n_r, n_c);
	Image resg(gear.n_rows, gear.n_cols);
	for (uint i = 0; i < n_r; i++)
		for (uint j = 0; j < n_c; j++){
			if (lbl_Img(i, j))
				res(i, j) = make_tuple(255, 255, 255);
			else
				res(i, j) = make_tuple(0, 0, 0);
			//if (i == static_cast<unsigned>(posg[num_dot].first) && j == static_cast<unsigned>(posg[num_dot].second))
			//	res(i, j) = make_tuple(255, 0, 0);
		}
	for (uint i = 0; i < gear.n_rows; i++)
		for (uint j = 0; j < gear.n_cols; j++){
			if (binImgGear(i, j))
				resg(i, j) = make_tuple(255, 255, 255);
			else
				resg(i, j) = make_tuple(0, 0, 0);
			//if (i == static_cast<unsigned>(posg[1].first) && j == static_cast<unsigned>(posg[1].second))
			//	res(i, j) = make_tuple(255, 0, 0);
		}
	return make_tuple(result_idx, object_array, res, resg);
}

tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in)
{
    // Base: return array of found objects and index of the correct gear
    // Bonus: return additional parameters of gears
	Matrix< uint > bn(in.n_cols, in.n_rows); //обычная матрица беззнаковых чисел

    auto object_array = vector<shared_ptr<IObject>>();
    int result_idx = 0;
    return make_tuple(result_idx, object_array, in.deep_copy());
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }
	/*
		загрузка изображения 
		создание вектора умных указателей на объект типа кортеж int int
		создание изображения, индекс результата
		создание кортежа из рез.индекс, вектора, изображ.
		сохранение
	*/
    try {
        Image src_image = load_image(argv[1]);
	Image src_gear = load_image(argv[4]);
        ofstream fout(argv[3]);
        vector<shared_ptr<IObject>> object_array; //создает умные указатели на объект типа кортеж int int
        Image dst_image, gear_image;
        int result_idx;
        //tie(result_idx, object_array, dst_image) = repair_mechanism(src_image);
	tie(result_idx, object_array, dst_image, gear_image) = repaire(src_image, src_gear);
        save_image(dst_image, argv[2]);
	save_image(gear_image, "out_gear.bmp");
        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
