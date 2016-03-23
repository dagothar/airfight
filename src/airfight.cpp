#include <allegro.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

//#define TEST1
#define TEST2
#define WRAP
#define FPS

const int xres = 800;
const int yres = 600;

const int mapx = 1600;
const int mapy = 1200;

int timer, etimer;

PALETTE main_pal;
BITMAP *bufor, *map, *port, *l_sky, *sky, *dmgmap, *fight, *ready, *splash;
BITMAP *speedm, *explosion, *ddoll, *bddoll;
BITMAP *shot, *cannon, *rocket, *bomb, *smoke;

std::string c1, c2;

enum buttons { NOWA, OPCJE, ZAKONCZ, WROC, W_1, LISTA1, W_2, LISTA2, T_PVPTIME, O_PVPTIME, T_FPS, O_FPS, T_DMGMUL, O_DMGMUL, C_FPS, C_WRAP, C_TEST2, B_OOK, B_DEFAULT, BTERM, BEND };
DIALOG m[BEND];

/* stałe */
const int max_players = 2;
const int portw = 390;
const int porth = 540;
const double PI = 3.141592;
const double dt = 0.01;
double at = 0.01;
const double scale = 10.0;
const double fuel_energy = 10000.0;
const double air_drag = 2.5;
const double THR_FORWARD = 1.0;
const double STR_LEFT = -1.0;
const double STR_RIGHT = 1.0;
const double max_leak = 1.0;
const double max_jam_chance = 0.75;
double no_pvp_time = 5.0;
double damage_mul = 1.0;
const double c_boost = 1.5;

std::string defaultdata("data/spitfire.dat");

/* flagi główne */
int xsize, ysize;
char f_menu = 1;
char quit = 0;
char game = 0;
double col_r = 0.0;
char wpn1 = 0, wpn2 = 0;
char bst1 = 0, bst2 = 0;
int p1score = 0, p2score = 0;
char engage = 0;
int frame = 0;
int f_timer;
int fps = 0;
char show_fps = 0;
char map_wraps = 0;
char test2 = 0;
volatile char go_pass = 0;
volatile char display_pass = 0;

/* samolot */
enum Shots { GUN, CANNON, ROCKET, BOMB };
struct Gun {
	int dx, dy;
	double dfi;
};

struct Weapon {
	std::string name;
	int type;
	int r0, r;
	double s_power, s_delay;
	double b_speed, b_life;
	int n_guns;
	struct Gun *guns_list;
};

enum Dzones { HULL, ENGINE, CREW, TANK, WEAPONS, SLEFT, SRIGHT, DEND, ALL};
enum Player_IDs { PLAYER1, PLAYER2 };

class Aircraft {
	private:
		char player;
		std::string name;
		char aloft, alive, crew_alive;
		double x, y, fi;
		int px, py;
		double v;
		double m;
		double f0, f;
		double eng_power0, idle_flow, leak;
		double eng_power, booster;
		double l_steer0, r_steer0, min_v, max_v;
		double l_steer, r_steer;
		double thr, str;
		char boost;
		int s_timer;
		PALETTE pal;
		BITMAP *sprite, *dmgsprite;
		double size;
		double hps0[DEND], hps[DEND];
		double smokep;

		int n_weapons, c_weapon;
		double jam_chance;
		struct Weapon *weapons_list;

		void init();
		void paint_plane();
	public:
		Aircraft(char _p, double _x, double _y, double _fi, std::string filename);
		void drawdoll(int ix, int iy);
		void info(int ix, int iy);
		void display();
		void ex_display();
		void go();
		void damage(double dmg, int area);
		void die();

		void set_thr(double set) { if(crew_alive) thr = set; }
		void set_str(double set) { if(crew_alive) str = set; }
		void wpn_change() { if(crew_alive) ++c_weapon; if(c_weapon >= n_weapons) c_weapon = 0; }
		void fire();
		void toogle_boost() { boost = !boost; }

		double get_x() { return x; }
		double get_y() { return y; }
		int get_px() { return px; }
		int get_py() { return py; }
		double get_size() { return size; }
		int fired, hit;
} *p1 = NULL, *p2 = NULL, *p[max_players];

/* pociski */
class Shot {
	private:
		double x, y, fi, v;
		double last_e, ei;
		int px, py;
		double power, b_life;
		int self_timer;
		double S, range;
		char owner;
	public:
		Shot(int _t, double _x, double _y, double _v, double _fi, double _power, double _life, char _owner);

		void display();
		void go();
		int check_hit();
		void die();

		int type;
		char active;
		Shot *next;
} *s = NULL;

/* kolory */
enum Colors { WHITE, BLACK, RED, GREEN, BLUE, YELLOW, TRACE1, TRACE2, CHULL, CENGINE, CCREW, CTANK, CWEAPONS, CSLEFT, CSRIGHT, TEXT, HUD, _END_ };
int color[_END_];

void make_colors()
{
	color[WHITE] = makecol(255, 255, 255);
	color[BLACK] = makecol(0, 0, 0);
	color[RED] = makecol(255, 0, 0);
	color[BLUE] = makecol(0, 0, 255);
	color[GREEN] = makecol(0, 150, 0);
	color[YELLOW] = makecol(255, 255, 0);
	color[TRACE1] = makecol(240, 240, 240);
	color[TRACE2] = makecol(150, 150, 255);

	color[CENGINE] = makecol(255, 0, 0);
	color[CCREW] = makecol(255, 128, 0);
	color[CTANK] = makecol(0, 255, 0);
	color[CHULL] = makecol(0, 128, 0);
	color[CWEAPONS] = makecol(255, 255, 0);
	color[CSLEFT] = makecol(128, 128, 0);
	color[CSRIGHT] = makecol(128, 0, 0);

	color[TEXT] = makecol(150, 150, 150);
	color[HUD] = makecol(50, 50, 50);
}

bool collide_sprite()
{
}

Aircraft::Aircraft(char _p, double _x, double _y, double _fi, std::string filename)
{
	player = _p; x = _x; y = _y; fi = _fi;

	/* wczytaj dane z pliku */
	std::ifstream data(filename.c_str(), std::ios::in);
	if(!data) data.open(defaultdata.c_str(), std::ios::in);

	/* wczytaj sprite */
	std::string spritefile;
	getline(data, spritefile);
	sprite = load_bitmap(spritefile.c_str(), pal);

	std::string dmgfile;
	getline(data, dmgfile);
	dmgsprite = load_bitmap(dmgfile.c_str(), pal);

	/* wczytaj parametry */
	std::string raw_para;
	getline(data, raw_para);
	std::stringstream para(raw_para);
	para >> name >> m >> f0 >> eng_power0 >> booster >> l_steer0 >> r_steer0 >> min_v >> max_v;

	/* wczytaj uzbrojenie */
	getline(data, raw_para);
	para.clear();
	para.str(raw_para);
	para >> n_weapons;
	weapons_list = new struct Weapon[n_weapons]; // ilosc broni

	int i, j;
	for(i = 0; i < n_weapons; ++i) {
		getline(data, raw_para);
		para.clear();
		para.str(raw_para);
		para >> weapons_list[i].name >> weapons_list[i].type >> weapons_list[i].r0 >> weapons_list[i].s_power >> weapons_list[i].s_delay >> weapons_list[i].b_speed >> weapons_list[i].b_life >> weapons_list[i].n_guns;
		weapons_list[i].r = weapons_list[i].r0;
		weapons_list[i].guns_list = new struct Gun[weapons_list[i].n_guns];
		getline(data, raw_para);
		para.clear();
		para.str(raw_para);
		for(j = 0; j < weapons_list[i].n_guns; ++j)
			para >> weapons_list[i].guns_list[j].dx >> weapons_list[i].guns_list[j].dy >> weapons_list[i].guns_list[j].dfi;
	}

	/* opancerzenie itp. */
	getline(data, raw_para);
	para.clear();
	para.str(raw_para);
	para >> size;
	for(i = 0; i < DEND; para >> hps0[i++]);

	data.close();

	init();
}


void Aircraft::paint_plane()
{
	int i, j;
	int c, n;
	int m = player;
	if(player == PLAYER1) c = color[RED];
	else if(player == PLAYER2) c = color[BLUE];
	for(i = 0; i < sprite->w; ++i)
	for(j = 0; j < sprite->h; ++j) {
		if(((long *)sprite->line[j])[i] == color[YELLOW]) ((long *)sprite->line[j])[i] = c;
		if(((long *)dmgsprite->line[j])[i] != 0) {
			n = ((long *)dmgsprite->line[j])[i];
			n = makecol(getr32(n), getg32(n), m);
			((long *)dmgsprite->line[j])[i] = n;
		}
	}
}

void Aircraft::init()
{
	aloft = 0; alive = 1; crew_alive = 1;
	v = 0.0;
	f = f0;
	thr = str = 0.0;
	eng_power = eng_power0;
	l_steer = l_steer0;
	r_steer = r_steer0;
	idle_flow = 0.25 * eng_power * dt;
	leak = 0.0;
	int i;
	for(i = 0; i < DEND; hps[i] = hps0[i++]);
	c_weapon = 0;
	boost = 0;
	jam_chance = 0.0;
	s_timer = clock();
	fired = hit = 0;
	smokep = 0.0;
	paint_plane();
}

void Aircraft::drawdoll(int ix, int iy)
{
	int i, j, k;
	int c[7];
	for(i = 0; i < DEND; ++i) c[i] = makecol((1.0 - hps[i] / hps0[i]) * 255, 0, 0);

	blit(ddoll, bddoll, 0, 0, 0, 0, ddoll->w, ddoll->h);

	for(i = 0; i < bddoll->w; ++i)
	for(j = 0; j < bddoll->h; ++j) {
		k = ((long *)bddoll->line[j])[i];
		if(k == color[CENGINE]) ((long *)bddoll->line[j])[i] = c[ENGINE];
		else if(k == color[CCREW]) ((long *)bddoll->line[j])[i] = c[CREW];
		else if(k == color[CTANK]) ((long *)bddoll->line[j])[i] = c[TANK];
		else if(k == color[CWEAPONS]) ((long *)bddoll->line[j])[i] = c[WEAPONS];
		else if(k == color[CSRIGHT]) ((long *)bddoll->line[j])[i] = c[SRIGHT];
		else if(k == color[CSLEFT]) ((long *)bddoll->line[j])[i] = c[SLEFT];
		else  if(k == color[CHULL]) ((long *)bddoll->line[j])[i] = c[HULL];
	}
	masked_blit(bddoll, bufor, 0, 0, ix, iy, ddoll->w, ddoll->h);
}

void Aircraft::info(int ix, int iy)
{
	textprintf_ex(bufor, font, ix, iy, color[TEXT], -1,
		"F");
	rect(bufor, ix+9, iy-1, ix+61, iy+11, color[TEXT]);
	rectfill(bufor, ix+10, iy, ix+(f/f0)*50+10, iy+10, color[RED]);

	double a = 1.0 * weapons_list[c_weapon].r / weapons_list[c_weapon].r0;
	textprintf_ex(bufor, font, ix+65, iy, color[TEXT], -1,
		"A");
	rect(bufor, ix+74, iy-1, ix+126, iy+11, color[TEXT]);
	rectfill(bufor, ix+75, iy, ix+a*50+75, iy+10, color[YELLOW]);

	textprintf_ex(bufor, font, ix+130, iy, color[TEXT], -1,
		"H");
	rect(bufor, ix+139, iy-1, ix+191, iy+11, color[TEXT]);
	rectfill(bufor, ix+140, iy, ix+(hps[HULL]/hps0[HULL])*50+140, iy+10, color[GREEN]);

	double alfa = -PI/2 + (v/max_v)*PI;
	masked_blit(speedm, bufor, 0, 0, ix+195, iy-11, speedm->w, speedm->h);
	line(bufor, ix+211, iy+5, ix+211+10.0*sin(alfa), iy+5-10.0*cos(alfa), color[BLUE]);

	textprintf_ex(bufor, font, ix+230, iy, color[TEXT], -1,
		"%s", weapons_list[c_weapon].name.c_str());

	//textprintf_ex(bufor, font, ix+300, iy, color[WHITE], -1,
	//	"Score: %d Accuracy: %.0f%%", player == 0 ? p1score : p2score, 100.0*hit/fired);
	textprintf_ex(bufor, font, ix, iy+15, color[TEXT], -1,
		"%s Booster: %.2f", name.c_str(), booster);
	textprintf_ex(bufor, font, ix+230, iy+15, color[TEXT], -1,
		"Score: %d", player == 0 ? p1score : p2score);

	textprintf_ex(bufor, font, ix+230, iy-15, color[TEXT], -1,
		"%2.0f %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f", hps[0], hps[1], hps[2], hps[3], hps[4], hps[5], hps[6]);

	drawdoll(ix+325, iy-5);
}

void Aircraft::display()
{
	if(alive) {
		int clr;
		if(player == PLAYER1) clr = color[TRACE1];
		else if(player == PLAYER2) clr = color[TRACE2];
		putpixel(sky, px, py, clr);
		if((1.0*rand()/RAND_MAX) < smokep)
			rotate_sprite(sky, smoke, px-smoke->w/2, py-smoke->h/2, fixmul(ftofix(fi), radtofix_r));
		rotate_sprite(bufor, sprite, px-sprite->w/2, py-sprite->h/2, fixmul(ftofix(fi), radtofix_r));
		rotate_sprite(dmgmap, dmgsprite, px-sprite->w/2, py-sprite->h/2, fixmul(ftofix(fi), radtofix_r));
	}
}

void Aircraft::ex_display()
{
	if(alive) {
		int clr;
		if(player == PLAYER1) clr = color[TRACE1];
		else if(player == PLAYER2) clr = color[TRACE2];
		putpixel(sky, px, py, clr);
		if((1.0*rand()/RAND_MAX) < smokep)
			rotate_sprite(sky, smoke, px-smoke->w/2, py-smoke->h/2, fixmul(ftofix(fi), radtofix_r));
		rotate_sprite(map, sprite, px-sprite->w/2, py-sprite->h/2, fixmul(ftofix(fi), radtofix_r));
		//rotate_sprite(dmgmap, dmgsprite, px-sprite->w/2, py-sprite->h/2, fixmul(ftofix(fi), //radtofix_r));
	}
}

void Aircraft::go()
{
	if(alive && hps[HULL] <= 0.0) die();
	if(!alive) return;
	if(booster <= 0.0) booster = 0.0, boost = 0;
	if(boost) booster -= dt;
	double flow = (eng_power * (boost ? c_boost : 1.0)) * abs(thr) * dt;
	f -= idle_flow + leak;
	if(f < 0.0) f = 0.0;
	flow = flow > f ? f : flow;
	f -= flow;
	if(f < 0.0) f = 0.0;
	double a = 0.0;
	if(v < max_v * (boost ? c_boost : 1.0)) a = (thr * flow * fuel_energy);
	a -= air_drag * v;
	v += a / m;
	if(v > min_v) aloft = 1;
	if(aloft && v < min_v) damage(1000.0, HULL);

	a = (v/(max_v-min_v));
	if(a < 0.0) a = 0.0;
	fi += ((str > 0.0) ? r_steer : l_steer) * str * a * dt;

	if(fi > PI) fi = fi - 2 * PI;
	if(fi < -PI) fi = 2 * PI + fi;
	x += v * sin(fi) * dt;
	y += v * cos(fi) * dt;
	px = xsize/2 + x * scale;
	py = ysize/2 - y * scale;

	#ifdef WRAP
	if(map_wraps) {
		if(x > (xsize/2) / scale) x -= xsize / scale;
		if(y > (ysize/2) / scale) y -= ysize / scale;
		if(x < -(xsize/2) / scale) x += xsize / scale;
		if(y < -(ysize/2) / scale) y += ysize / scale;
	} else {
		if(x > (xsize/2) / scale) damage(1000.0, HULL);
		if(y > (ysize/2) / scale) damage(1000.0, HULL);
		if(x < -(xsize/2) / scale) damage(1000.0, HULL);
		if(y < -(ysize/2) / scale) damage(1000.0, HULL);
	}
	#endif

	if(!crew_alive) {
		str += 0.1 * ((2.0 * rand() / RAND_MAX) - 1.0);
		if(str < STR_LEFT) str = STR_LEFT;
		if(str > STR_RIGHT) str = STR_RIGHT;
		thr += 0.1 * ((2.0 * rand() / RAND_MAX) - 1.0);
		if(thr < 0.0) thr = 0.0;
		if(thr > THR_FORWARD) thr = THR_FORWARD;
	} else {
		#ifdef TEST2
		if((player == PLAYER2) && crew_alive && test2) {
			str += 0.1 * ((2.0 * rand() / RAND_MAX) - 1.0);
			if(str < STR_LEFT) str = STR_LEFT;
			if(str > STR_RIGHT) str = STR_RIGHT;
			thr = THR_FORWARD;
		} else {
		#endif
			str = 0.0;
			thr = 0.0;
		#ifdef TEST2
		}
		#endif
	}
	rotate_sprite(dmgmap, dmgsprite, px-sprite->w/2, py-sprite->h/2, fixmul(ftofix(fi), radtofix_r));
}

void Aircraft::fire()
{
	if((clock() - s_timer > weapons_list[c_weapon].s_delay * CLOCKS_PER_SEC) && weapons_list[c_weapon].r > 0 && alive && aloft && crew_alive && engage) {
		int n;
		double sx, sy, dfi, dx, dy;
		for(n = 0; n < weapons_list[c_weapon].n_guns && weapons_list[c_weapon].r > 0 && (1.0 * rand() / RAND_MAX) > jam_chance; ++n) {
			sx = x-((sprite->w/2-weapons_list[c_weapon].guns_list[n].dx)*cos(fi)-(sprite->h/2-weapons_list[c_weapon].guns_list[n].dy)*sin(fi))/scale;
			sy = y+((sprite->h/2-weapons_list[c_weapon].guns_list[n].dy)*cos(fi)+(sprite->w/2-weapons_list[c_weapon].guns_list[n].dx)*sin(fi))/scale;
			dfi = weapons_list[c_weapon].guns_list[n].dfi;
			if(dfi > PI) {
				dx = p1->get_x() - p2->get_x();
				dy = p1->get_y() - p2->get_y();
				dfi = atan2(dx, dy) + PI;
				if(player == PLAYER2) dfi -= PI;
			}
			else dfi += fi;
			Shot *new_shot = new Shot(weapons_list[c_weapon].type, sx, sy, v+weapons_list[c_weapon].b_speed, dfi, weapons_list[c_weapon].s_power, weapons_list[c_weapon].b_life, player);
			new_shot->next = s;
			s = new_shot;
			++fired;
			--weapons_list[c_weapon].r;
		}
		s_timer = clock();
	}
}

void Aircraft::damage(double dmg, int area)
{
 //WHITE, HULL, ENGINE, CREW, TANK, WEAPONS, SLEFT, SRIGHT
	int rm;
	dmg *= damage_mul;
	switch(area)
	{
		case ENGINE:
			hps[ENGINE] -= dmg;
			if(hps[ENGINE] <= 0.0) hps[ENGINE] = 0.0;
			eng_power = eng_power0 * hps[ENGINE] / hps0[ENGINE];
			break;
		case CREW:
			hps[CREW] -= dmg;
			//dmg *= 2;
			if(hps[CREW] <= 0.0) hps[CREW] = 0.0, crew_alive = 0;
			break;
		case TANK:
			hps[TANK] -= dmg;
			if(hps[TANK] <= 0.0) hps[TANK] = 0.0, die();
			leak = max_leak * (1.0 - hps[TANK] / hps0[TANK]);
			if(smokep < 0.1 * leak) smokep = 0.1 * leak / max_leak;
			break;
		case WEAPONS:
			hps[WEAPONS] -= dmg;
			if(hps[WEAPONS] <= 0.0) hps[WEAPONS] = 0.0;
			rm = dmg * weapons_list[c_weapon].r0 / hps0[WEAPONS];
			if(rm < 1) rm = 1;
			weapons_list[c_weapon].r -= rm;
			if(weapons_list[c_weapon].r < 0) weapons_list[c_weapon].r = 0;
			jam_chance = max_jam_chance * (1.0 - hps[WEAPONS] / hps0[WEAPONS]);
			break;
		case SLEFT:
			hps[SLEFT] -= dmg;
			if(hps[SLEFT] <= 0.0) hps[SLEFT] = 0.0;
			l_steer = l_steer0 * (hps[SLEFT] / hps0[SLEFT]);
			break;
		case SRIGHT:
			hps[SRIGHT] -= dmg;
			if(hps[SRIGHT] <= 0.0) hps[SRIGHT] = 0.0;
			r_steer = r_steer0 * (hps[SRIGHT] / hps0[SRIGHT]);
			break;
		case HULL:
			hps[HULL] -= 0.5 * dmg;
			break;
	}
	hps[HULL] -= 0.5 * dmg;
	if(hps[HULL] <= 0.0) hps[HULL] = 0.0;
}

void Aircraft::die()
{
	v = 0.0;
	eng_power = 0.0;
	l_steer = r_steer = 0.0;
	f = 0.0;
	col_r = 0.0;
	size = 0.0;
	alive = 0;

	if(player == PLAYER1) ++p2score;
	if(player == PLAYER2) ++p1score;

	rotate_sprite(sky, explosion, px-explosion->w/2, py-explosion->h/2, fixmul(ftofix(fi), radtofix_r));
}

Shot::Shot(int _t, double _x, double _y, double _v, double _fi, double _power, double _life, char _owner)
{
	type = _t;
	x = _x; y = _y; fi = _fi; power = _power;
	b_life = _life;
	v = _v;
	px = xsize/2 + x * scale; py = ysize/2 - scale * y;
	self_timer = clock();
	next = NULL;
	active = 1;
	last_e = ei = 0.0;
	owner = _owner;
}

void Shot::display()
{
	BITMAP *spr;
	switch(type) {
		case GUN: spr = shot; break;
		case CANNON: spr = cannon; break;
		case ROCKET: spr = rocket; putpixel(sky, px, py, color[TRACE1]); break;
		case BOMB: spr = bomb; break;
		default: spr = shot; break;
	}
	rotate_sprite(map, spr, px - spr->w/2, py - spr->h/2, fixmul(ftofix(fi), radtofix_r));
}

void Shot::go()
{
	if(type == BOMB) v -= air_drag * 0.1 * v * dt;

	double dmg = (type == BOMB ? 0.0 : power);
	int hit = check_hit();
	if(hit)
		if((hit & 0x08) && owner != PLAYER2) {
			p2->damage(dmg, hit & 0x07);
			die();
			if(owner != -1) ++p1->hit;
		}
		else if(!(hit & 0x08) && owner != PLAYER1) {
			p1->damage(dmg, hit & 0x07);
			die();
			if(owner != -1) ++p2->hit;
		}

	x += v * sin(fi) * dt;
	y += v * cos(fi) * dt;
	px = xsize/2 + x * scale;
	py = ysize/2 - y * scale;
	#ifdef WRAP
	if(map_wraps) {
		if(x > (xsize/2) / scale) x -= xsize / scale;
		if(y > (ysize/2) / scale) y -= ysize / scale;
		if(x < -(xsize/2) / scale) x += xsize / scale;
		if(y < -(ysize/2) / scale) y += ysize / scale;
	} else {
		if(x > (xsize/2) / scale) die();
		if(y > (ysize/2) / scale) die();
		if(x < -(xsize/2) / scale) die();
		if(y < -(ysize/2) / scale) die();
	}
	#endif

	double dx = x, dy = y;
	if(owner == PLAYER1) {
		dx -= p2->get_x();
		dy -= p2->get_y();
	} else if(owner == PLAYER2) {
		dx -= p1->get_x();
		dy -= p1->get_y();
	}

	if(type == ROCKET) {
		double cfi = atan2(-dx, -dy), efi, e1, e2, dfi;
		e1 = fi - cfi;
		e2 = PI-cfi-(-PI-fi);
		efi = abs(e1) < abs(e2) ? e1 : e2;
		cfi = v / 10 * dt;
		dfi = cfi * (efi + 0.1 * (efi - last_e) / dt + 0.05 * ei);
		if(dfi > cfi) dfi = cfi;
		if(dfi < -cfi) dfi = -cfi;
		fi -= dfi;

		if(fi > PI) fi = fi - 2 * PI;
		if(fi < -PI) fi = 2 * PI + fi;
		last_e = efi;
		ei += efi * dt;
	}

	if(clock() - self_timer > b_life * CLOCKS_PER_SEC) die();
}

int Shot::check_hit()
{
	int c, r, o = 0;
	if(is_inside_bitmap(dmgmap, px, py, 0)) {
		c = ((long *)dmgmap->line[py])[px];
		o = getb32(c);
		if(c != 0 && o != owner) {
			c = makecol(getr32(c), getg32(c), 0);
			if(c == color[CENGINE]) r = ENGINE;
			else if(c ==color[CCREW]) r = CREW;
			else if(c == color[CTANK]) r = TANK;
			else if(c == color[CWEAPONS]) r = WEAPONS;
			else if(c == color[CSRIGHT]) r = SRIGHT;
			else if(c == color[CSLEFT]) r = SLEFT;
			else if(c == color[CHULL]) r = HULL;
			r ^= o << 3;
			r ^= 0xf0;
			return r;
		}
	}
	return 0;
}

void Shot::die()
{
	active = 0;
	double f;
	switch(type) {
		case ROCKET: rotate_sprite(sky, explosion, px-explosion->w/2, py-explosion->h/2, fixmul(ftofix(fi), radtofix_r)); return;
		case CANNON:
			rotate_sprite(sky, smoke, px-smoke->w/2, py-smoke->h/2, fixmul(ftofix(fi), radtofix_r));
			for(f = -PI; f < PI; f += PI / (rand() % 3 + 3)) {
				Shot *new_shot = new Shot(GUN, x, y, v+75.0, f, power/2, 0.01, -1);
				new_shot->next = s;
				s = new_shot;
			}
			return;
		case BOMB:
			rotate_sprite(sky, explosion, px-explosion->w/2, py-explosion->h/2, fixmul(ftofix(fi), radtofix_r));
			for(f = -PI; f < PI; f += PI / (rand() % 10 + 10)) {
				Shot *new_shot = new Shot(GUN, x, y, v+50.0, f, power, 0.15, -1);
				new_shot->next = s;
				s = new_shot;
			}
			return;
		default: return;
	}
}

void purge_bullets()
{
	Shot *t_s = s;
	while(t_s) {
		t_s->active = 0;
		t_s = t_s->next;
	}
}

void place_aircraft(std::string c1, std::string c2)
{
	engage = 0;

	purge_bullets();

	if(p1) delete p1;
	if(p2) delete p2;

	p1 = new Aircraft(PLAYER1, -5.0, -mapy / (2*scale) + 10.0, 0.0, "data/"+c1+".dat");
	p2 = new Aircraft(PLAYER2, 5.0, mapy / (2*scale) - 10.0, PI, "data/"+c2+".dat");

	col_r = (p1->get_size() > p2->get_size()) ? p1->get_size() : p2->get_size();

	etimer = clock();
}

void go_timer()
{
	go_pass = 1;
}
END_OF_FUNCTION(go_timer);

void display_timer()
{
	display_pass = 1;
}
END_OF_FUNCTION(display_timer);

void init()
{
	allegro_init();
	install_timer();
	install_keyboard();
	install_mouse();
	srand(time(NULL));

	set_color_depth(32);
	set_gfx_mode(GFX_AUTODETECT_WINDOWED, xres, yres, 0, 0);

	set_window_title("Airfight");

	make_colors();

	bufor = create_bitmap(xres, yres);
	clear_bitmap(bufor);
	sky = create_bitmap(mapx, mapy);
	clear_bitmap(sky);
	l_sky = load_bitmap("data/sky.bmp", main_pal);
	stretch_blit(l_sky, sky, 0, 0, l_sky->w, l_sky->h, 0, 0, mapx, mapy);

	xsize = mapx;
	ysize = mapy;
	map = create_bitmap(xsize, ysize);
	clear_bitmap(map);
	dmgmap = create_bitmap(xsize, ysize);
	clear_bitmap(dmgmap);

	splash = load_bitmap("data/menu.bmp", main_pal);
	speedm = load_bitmap("data/speed.bmp", main_pal);
	ddoll = load_bitmap("data/ddoll.bmp", main_pal);
	bddoll = create_bitmap(ddoll->w, ddoll->h);
	fight = load_bitmap("data/fight.bmp", main_pal);
	ready = load_bitmap("data/ready.bmp", main_pal);
	explosion = load_bitmap("data/explosion.bmp", main_pal);
	shot = load_bitmap("data/shot.bmp", main_pal);
	cannon = load_bitmap("data/cannon.bmp", main_pal);
	rocket = load_bitmap("data/rocket.bmp", main_pal);
	bomb = load_bitmap("data/bomb.bmp", main_pal);
	smoke = load_bitmap("data/smoke.bmp", main_pal);

	install_int_ex(go_timer, MSEC_TO_TIMER(dt * 1000));
	LOCK_VARIABLE(go_pass);
	LOCK_FUNCTION(go_timer);
	install_int_ex(display_timer, MSEC_TO_TIMER(at * 1000));
	LOCK_VARIABLE(display_pass);
	LOCK_FUNCTION(display_timer);
	f_timer = clock();
}

void check_controls()
{
	poll_keyboard();
	/* klawisze główne */
	if(key[KEY_ESC]) { f_menu = 1; return; }
	if(key[KEY_P]) { place_aircraft(c1, c2); return; }

	/* klawisze gracza #1 */
	if(key[KEY_UP]) p1->set_thr(THR_FORWARD);
	if(key[KEY_DOWN] && !bst1) bst1 = 1, p1->toogle_boost();
	if(!key[KEY_DOWN]) bst1 = 0;
	if(key[KEY_LEFT]) p1->set_str(STR_LEFT);
	if(key[KEY_RIGHT]) p1->set_str(STR_RIGHT);
	if(key[KEY_RCONTROL]) p1->fire();
	if(key[KEY_RSHIFT] && !wpn1) wpn1 = 1, p1->wpn_change();
	if(!key[KEY_RSHIFT]) wpn1 = 0;

	/* klawisze gracza #2 */
	if(key[KEY_X]) p2->set_thr(THR_FORWARD);
	if(key[KEY_A] && !bst2) bst2 = 1, p2->toogle_boost();
	if(!key[KEY_A]) bst2 = 0;
	if(key[KEY_Z]) p2->set_str(STR_LEFT);
	if(key[KEY_C]) p2->set_str(STR_RIGHT);
	if(key[KEY_LCONTROL]) p2->fire();
	if(key[KEY_LSHIFT] && !wpn2) wpn2 = 1, p2->wpn_change();
	if(!key[KEY_LSHIFT]) wpn2 = 0;

	#ifdef TEST1
	if(test1) {
		p2->set_thr(THR_FORWARD);
		p2->set_str(STR_RIGHT);
	}
	#endif

	clear_keybuf();
}

void display()
{
	clear_to_color(bufor, color[HUD]);

	int x1 = p1->get_px()-portw/2;
	int x2 = p2->get_px()-portw/2;
	int y1 = p1->get_py()-porth/2;
	int y2 = p2->get_py()-porth/2;

	blit(sky, map, x1, y1, x1, y1, portw, porth);
	blit(sky, map, x2, y2, x2, y2, portw, porth);

	p1->info(400, yres-30);
	p2->info(0, yres-30);
	#ifdef FPS
	if(show_fps) {
		textprintf_ex(bufor, font, 0, 0, color[TEXT], 0,
		"FPS: %d", fps);
		++frame;
		if(clock() - f_timer > CLOCKS_PER_SEC) {
			fps = frame;
			frame = 0;
			f_timer = clock();
		}
	}
	#endif

	p1->ex_display();
	p2->ex_display();

	Shot *t_s = s;
	while(t_s) {
		t_s->display();
		t_s = t_s->next;
	}

	rectfill(bufor, xres/2+5, 10, xres/2+5+portw, 10+porth, color[BLACK]);
	rectfill(bufor, 5, 10, 5+portw, 10+porth, color[BLACK]);
	blit(map, bufor, x1, y1, xres/2+5, 10, portw, porth);
	blit(map, bufor, x2, y2, 5, 10, portw, porth);
	if(!engage) masked_blit(ready, bufor, 0, 0, xres/2-ready->w/2, yres/2-ready->h/2, ready->w, ready->h);
	else if(clock() - etimer < no_pvp_time * CLOCKS_PER_SEC)
		masked_blit(fight, bufor, 0, 0, xres/2-fight->w/2, yres/2-fight->h/2, fight->w, fight->h);

	blit(bufor, screen, 0, 0, 0, 0, xres, yres);
}

void go()
{
	if((!engage) && (clock() - etimer > no_pvp_time * CLOCKS_PER_SEC)) {
		engage = 1;
		etimer = clock();
	}

	rectfill(dmgmap, p1->get_px()-25, p1->get_py()-25, p1->get_px()+25, p1->get_py()+25, 0);
	rectfill(dmgmap, p2->get_px()-25, p2->get_py()-25, p2->get_px()+25, p2->get_py()+25, 0);

	p1->go();
	p2->go();

	/* sprawdz, czy kolizja */
	double R, dx, dy;
		dx = p1->get_x() - p2->get_x();
		dy = p1->get_y() - p2->get_y();
		R = dx * dx + dy * dy;
		R = sqrt(R);
		if(R < col_r) p1->damage(1000.0, HULL), p2->damage(1000.0, HULL);


	/* pociski */
	Shot *t_s = s, *t_o = NULL, *t_next = NULL;;
	s = NULL;

	while(t_s) {
		if(t_s->active) {
			if(!s) s = t_s;
			t_s->go();
			t_o = t_s;
			t_s = t_s->next;
		} else {
			t_next = t_s->next;
			if(t_o) t_o->next = t_next;
			//t_s->die();
			delete t_s;
			t_s = t_next;
		}
	}
}

void new_game()
{
	stretch_blit(l_sky, sky, 0, 0, l_sky->w, l_sky->h, 0, 0, mapx, mapy);

	switch(m[LISTA1].d1) {
		case 0: c1 = "spitfire"; break;
		case 1: c1 = "bf110"; break;
		case 2: c1 = "b17"; break;
		case 3: c1 = "f16"; break;
		default: c1 = "spitfire";
	}

	switch(m[LISTA2].d1) {
		case 0: c2 = "spitfire"; break;
		case 1: c2 = "bf110"; break;
		case 2: c2 = "b17"; break;
		case 3: c2 = "f16"; break;
		default: c2 = "spitfire";
	}

	place_aircraft(c1, c2);

	game = 1;
}

int b_mbutton(int msg, DIALOG *d, int c)
{
	if(msg == MSG_CLICK || msg == MSG_KEY) {
		if(d == &m[ZAKONCZ]) {
			quit = 1;
			return D_EXIT;
		}
		if(d == &m[NOWA]) {
			new_game();
			f_menu = 0;
			return D_EXIT;
		}
		if(d == &m[WROC]) {
			f_menu = 0;
			return D_EXIT;
		}
		if(d == &m[OPCJE]) {
			m[NOWA].flags ^= D_DISABLED;
			m[WROC].flags ^= D_DISABLED;
			m[OPCJE].flags ^= D_DISABLED;
			m[ZAKONCZ].flags ^= D_DISABLED;
			m[LISTA1].flags ^= D_DISABLED;
			m[LISTA2].flags ^= D_DISABLED;

			m[O_PVPTIME].flags ^= D_HIDDEN;
			m[T_PVPTIME].flags ^= D_HIDDEN;
			m[O_DMGMUL].flags ^= D_HIDDEN;
			m[T_DMGMUL].flags ^= D_HIDDEN;
			m[O_FPS].flags ^= D_HIDDEN;
			m[T_FPS].flags ^= D_HIDDEN;
			m[C_FPS].flags ^= D_HIDDEN;
			m[C_WRAP].flags ^= D_HIDDEN;
			m[C_TEST2].flags ^= D_HIDDEN;
			m[B_OOK].flags ^= D_HIDDEN;
			m[B_DEFAULT].flags ^= D_HIDDEN;

			scare_mouse();
			rectfill(screen, 225, 75, xres-105, 430, color[HUD]);
			unscare_mouse();
			return D_REDRAW;
		}
		if(d == &m[B_OOK]) {
			m[NOWA].flags ^= D_DISABLED;
			m[WROC].flags ^= D_DISABLED;
			m[OPCJE].flags ^= D_DISABLED;
			m[ZAKONCZ].flags ^= D_DISABLED;
			m[LISTA1].flags ^= D_DISABLED;
			m[LISTA2].flags ^= D_DISABLED;

			m[O_PVPTIME].flags ^= D_HIDDEN;
			m[T_PVPTIME].flags ^= D_HIDDEN;
			m[O_DMGMUL].flags ^= D_HIDDEN;
			m[T_DMGMUL].flags ^= D_HIDDEN;
			m[O_FPS].flags ^= D_HIDDEN;
			m[T_FPS].flags ^= D_HIDDEN;
			m[C_FPS].flags ^= D_HIDDEN;
			m[C_WRAP].flags ^= D_HIDDEN;
			m[C_TEST2].flags ^= D_HIDDEN;
			m[B_OOK].flags ^= D_HIDDEN;
			m[B_DEFAULT].flags ^= D_HIDDEN;

			scare_mouse();
			blit(splash, screen, 0, 0, xres/2-splash->w/2, yres/2-splash->h/2, splash->w, splash->h);
			unscare_mouse();
			return D_REDRAW;
		}
		if(d == &m[B_DEFAULT]) {
			show_fps = 0;
			map_wraps = 0;
			test2 = 0;
			no_pvp_time = 5.0;
			damage_mul = 1.0;
			at = 0.01;

			m[O_PVPTIME].d2 = 5;
			m[O_DMGMUL].d2 = 5;
			m[O_FPS].d2 = 90;
			m[C_FPS].flags &= ~D_SELECTED;
			m[C_WRAP].flags &= ~D_SELECTED;
			m[C_TEST2].flags &= ~D_SELECTED;
			return D_REDRAW;
		}
		if(d == &m[C_FPS]) {
			show_fps = !show_fps;
			return d_check_proc(msg, d, c);
		}
		if(d == &m[C_WRAP]) {
			map_wraps = !map_wraps;
			return d_check_proc(msg, d, c);
		}
		if(d == &m[C_TEST2]) {
			test2 = !test2;
			return d_check_proc(msg, d, c);
		}
	}
	if(d == &m[C_FPS] || d == &m[C_WRAP] || d == &m[C_TEST2]) return d_check_proc(msg, d, c);
	return d_button_proc(msg, d, c);
}

char *b_list(int index, int *list_size)
{
	if(index < 0) {
		*list_size = 4;
		return NULL;
	} else {
		switch(index) {
			case 0: return "Spitfire";
			case 1: return "Bf110";
			case 2: return "B17";
			case 3: return "F16";
		}
	}
}

int s_pvptime(void *dp3, int d2)
{
	no_pvp_time = 1.0 * d2;
	return D_O_K;
}

int s_dmgmul(void *dp3, int d2)
{
	damage_mul = 0.5 + 0.1 * d2;
	return D_O_K;
}

int s_fps(void *dp3, int d2)
{
	at = 0.1 - d2 * 0.001;
	install_int_ex(display_timer, MSEC_TO_TIMER(at * 1000));
	return D_O_K;
}

void make1_menu()
{
}

void make_menu()
{
	int i = 0;

	i = NOWA;
	m[i].proc = b_mbutton;
	m[i].x = 100; m[i].y = 400;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'n';
	m[i].flags = D_O_K;
	m[i].dp = (void *)"&Nowa gra";

	i = OPCJE;
	m[i].proc = b_mbutton;
	m[i].x = 100; m[i].y = 450;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'o';
	m[i].flags = D_O_K;
	m[i].dp = (void *)"&Opcje";

	i = ZAKONCZ;
	m[i].proc = b_mbutton;
	m[i].x = 100; m[i].y = 500;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'w';
	m[i].flags = D_O_K;
	m[i].dp = (void *)"&Wyjdź";

	i = WROC;
	m[i].proc = b_mbutton;
	m[i].x = 100; m[i].y = 350;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'p';
	m[i].flags = D_O_K;
	m[i].dp = (void *)"&Powrót";

	i = W_1;
	m[i].proc = d_text_proc;
	m[i].x = 250; m[i].y = 440;
	m[i].fg = color[WHITE]; m[i].bg = -1;
	m[i].flags = D_O_K;
	m[i].dp = (void *)"Samolot #1";

	i = LISTA1;
	m[i].proc = d_list_proc;
	m[i].x = 250; m[i].y = 450;
	m[i].w = 100; m[i].h = 75;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 0;
	m[i].flags = D_O_K;
	m[i].dp = (void *)b_list;

	i = W_2;
	m[i].proc = d_text_proc;
	m[i].x = 400; m[i].y = 440;
	m[i].fg = color[WHITE]; m[i].bg = -1;
	m[i].flags = D_O_K;
	m[i].dp = (void *)"Samolot #2";

	i = LISTA2;
	m[i].proc = d_list_proc;
	m[i].x = 400; m[i].y = 450;
	m[i].w = 100; m[i].h = 75;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].flags = D_O_K;
	m[i].dp = (void *)b_list;

	i = T_PVPTIME;
	m[i].proc = d_text_proc;
	m[i].x = 250; m[i].y = 100;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[WHITE]; m[i].bg = -1;
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].dp = (void *)"Czas na rozgrzewkę (0 - 10 sek.)";

	i = O_PVPTIME;
	m[i].proc = d_slider_proc;
	m[i].x = 250; m[i].y = 110;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 10;
	m[i].d2 = 5;
	m[i].dp2 = (void *)s_pvptime;

	i = T_DMGMUL;
	m[i].proc = d_text_proc;
	m[i].x = 250; m[i].y = 140;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[WHITE]; m[i].bg = -1;
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].dp = (void *)"Współczynnik obrażeń (0.5 - 2.0)";

	i = O_DMGMUL;
	m[i].proc = d_slider_proc;
	m[i].x = 250; m[i].y = 150;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 15;
	m[i].d2 = 5;
	m[i].dp2 = (void *)s_dmgmul;

	i = T_FPS;
	m[i].proc = d_text_proc;
	m[i].x = 250; m[i].y = 180;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[WHITE]; m[i].bg = -1;
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].dp = (void *)"Odświeżanie grafiki (10 - 100)";

	i = O_FPS;
	m[i].proc = d_slider_proc;
	m[i].x = 250; m[i].y = 190;
	m[i].w = 200; m[i].h = 20;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 90;
	m[i].d2 = 90;
	m[i].dp2 = (void *)s_fps;

	i = C_FPS;
	m[i].proc = b_mbutton;
	m[i].x = 250; m[i].y = 220;
	m[i].w = 200; m[i].h = 10;
	m[i].fg = color[WHITE]; m[i].bg = color[HUD];
	m[i].key = 'f';
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 1;
	m[i].dp = (void *)"Pokaż &Fps";

	i = C_WRAP;
	m[i].proc = b_mbutton;
	m[i].x = 250; m[i].y = 240;
	m[i].w = 200; m[i].h = 10;
	m[i].fg = color[WHITE]; m[i].bg = color[HUD];
	m[i].key = 'z';
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 1;
	m[i].dp = (void *)"&Zawijanie mapy";

	i = C_TEST2;
	m[i].proc = b_mbutton;
	m[i].x = 250; m[i].y = 260;
	m[i].w = 200; m[i].h = 10;
	m[i].fg = color[WHITE]; m[i].bg = color[HUD];
	m[i].key = 't';
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].d1 = 1;
	m[i].dp = (void *)"&Test2";

	i = B_OOK;
	m[i].proc = b_mbutton;
	m[i].x = 250; m[i].y = 400;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'o';
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].dp = (void *)"&Ok";

	i = B_DEFAULT;
	m[i].proc = b_mbutton;
	m[i].x = 375; m[i].y = 400;
	m[i].w = 100; m[i].h = 25;
	m[i].fg = color[TEXT]; m[i].bg = color[HUD];
	m[i].key = 'd';
	m[i].flags = D_O_K | D_HIDDEN;
	m[i].dp = (void *)"&Domyślne";


	i = BTERM;
	m[i].proc = NULL;
}

void menu()
{
	if(game) m[WROC].flags &= ~D_HIDDEN; else m[WROC].flags |= D_HIDDEN;

	clear(screen);
	blit(splash, screen, 0, 0, xres/2-splash->w/2, yres/2-splash->h/2, splash->w, splash->h);
	do_dialog(m, -1);
}

int main(int argc, char *argv[])
{
	init();
	make_menu();

	while(true) {
		if(f_menu) menu();
		else {
			check_controls();

			if(go_pass) {
				go();
				go_pass = 0;
			}
			if(display_pass) {
				display();
				display_pass = 0;
			}
		}
		if(quit) break;
	};

	destroy_bitmap(bufor);
	destroy_bitmap(map);
	destroy_bitmap(l_sky);
	destroy_bitmap(sky);
	destroy_bitmap(dmgmap);
	destroy_bitmap(fight);
	destroy_bitmap(ready);
	destroy_bitmap(speedm);
	destroy_bitmap(explosion);
	destroy_bitmap(ddoll);
	destroy_bitmap(bddoll);
	destroy_bitmap(shot);
	destroy_bitmap(cannon);
	destroy_bitmap(rocket);
	destroy_bitmap(bomb);
	destroy_bitmap(smoke);
	destroy_bitmap(splash);

	allegro_exit();
	return 0;
}
END_OF_MAIN();
