/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */

#include "Advanced.h"
namespace LFHPrimitive{

extern map< void* , pair< uint32_t, uint32_t >  > man_in_the_middle;
myHashmap< uint32_t , pair< void*, uint32_t> > AliasBank;
myHashmap< const void* , uint32_t > AliasOf;
RunningStatistics< > lfhst_stats;


OrthoRef3 OrthoRef3::operator*(const OrthoRef3& other) const{ OrthoRef3 fout;
    switch(((r >> 3) & 7) | (other.r & 56)){
        case 0: case 1: case 2: case 3: case 4: case 5:
            fout.r = r ^ other.r; // simplest of all, inherit directed and flip mismatching directions
        break;
        case 9:
        case 8:fout.r = ((r&57)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 16:fout.r = (((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 27:
        case 24:fout.r = ((r&60)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 32:fout.r = (((r&3)<<1)|((r&4) >> 2)) ^ other.r; break;

        case 45:
        case 40:fout.r = ((r&58)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 18:fout.r = 48 ^ ((((r&6)>>1)|((r&1) << 2)) ^ other.r);break;
        case 36:fout.r = 48 ^ ((((r&3)<<1)|((r&4) >> 2)) ^ other.r);break;

        case 10:fout.r = 16^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 11:fout.r = 24^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 12:fout.r = 32^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 13:fout.r = 40^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;

        case 26:fout.r = 48^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 25:fout.r = 56^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 28:fout.r = 16^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 29:fout.r = 8^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;

        case 42:fout.r = 32^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 43:fout.r = 8^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 44:fout.r = 48^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 41:fout.r = 56^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;

        case 17:fout.r = 56^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 19:fout.r = 24^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 20:fout.r = (((r&6)>>1)|((r&1) << 2)) ^ (other.r & 7);break;
        case 21:fout.r = 8^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;

        case 33:fout.r = 56^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
        case 34:fout.r = ((((r&3)<<1)|((r&4) >> 2)) ^ (other.r& 7));break;
        case 35:fout.r = 8^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
        case 37:fout.r = 40^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
    }
    return fout;
}
OrthoRef3& OrthoRef3::operator*=(const OrthoRef3& other){ OrthoRef3 fout;
    switch(((r >> 3) & 7) | (other.r & 56)){ // 6x6 cases!
        case 0: case 1: case 2: case 3: case 4: case 5:
            r ^= other.r; // simplest of all, inherit directed and flip mismatching directions
        break;
        case 9:
        case 8:r = ((r&57)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 16:r = (((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 27:
        case 24:r = ((r&60)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 32:r = (((r&3)<<1)|((r&4) >> 2)) ^ other.r; break;
        case 45:
        case 40:r = ((r&58)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;

        case 18:r = 48 ^ ((((r&6)>>1)|((r&1) << 2)) ^ other.r);break;

        case 36:r = 48 ^ ((((r&3)<<1)|((r&4) >> 2)) ^ other.r);break;

        case 10:r = 16^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 11:r = 24^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 12:r = 32^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;
        case 13:r = 40^((r&1)|((r&2) << 1)| ((r&4) >> 1)) ^ other.r;break;

        case 26:r = 48^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 25:r = 56^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 28:r = 16^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;
        case 29:r = 8^((r&4)|((r&1) << 1)| ((r&2) >> 1)) ^ other.r;break;

        case 42:r = 32^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 43:r = 8^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 44:r = 48^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;
        case 41:r = 56^((r&2)|((r&1) << 2)| ((r&4) >> 2)) ^ other.r;break;

        case 17:r = 56^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 19:r = 24^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;
        case 20:r = (((r&6)>>1)|((r&1) << 2)) ^ (other.r & 7);break;
        case 21:r = 8^(((r&6)>>1)|((r&1) << 2)) ^ other.r;break;

        case 33:r = 56^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
        case 34:r = ((((r&3)<<1)|((r&4) >> 2)) ^ (other.r& 7));break;
        case 35:r = 8^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
        case 37:r = 40^(((r&3)<<1)|((r&4) >> 2)) ^ other.r;break;
    }
    return *this;
}
OrthoRef3& OrthoRef3::operator/=(const OrthoRef3& other){return (*this) *= other.mkInverse();}
OrthoRef3 OrthoRef3::operator/(const OrthoRef3& other)const{return (*this) * other.mkInverse();}
OrthoRef3 OrthoRef3::mkInverse() const{
    switch(r >> 3){
        case 1: return OrthoRef3((r & 0x39) | ((r >> 1) & 2) | ((r << 1) & 4));
        case 2: return OrthoRef3(((r >> 2) & 1) | ((r << 1) & 0x26));
        case 3: return OrthoRef3((r & 0x3C) | ((r >> 1) & 1) | ((r << 1) & 2));
        case 4: return OrthoRef3(((r >> 1) & 0x13) | ((r << 2) & 4));
        case 5: return OrthoRef3((r & 0x3A) | ((r >> 2) & 1) | ((r << 2) & 4));
        default: return OrthoRef3(r); // case 0 and illegal
    }
}
OrthoRef3& OrthoRef3::toinverse(){
    switch(r >> 3){
        case 1: r = (r & 0x39) | ((r >> 1) & 2) | ((r << 1) & 4); break;
        case 2: r = ((r >> 2) & 1) | ((r << 1) & 0x26); break;
        case 3: r = (r & 0x3C) | ((r >> 1) & 1) | ((r << 1) & 2); break;
        case 4: r = ((r >> 1) & 0x13) | ((r << 2) & 4); break;
        case 5: r = (r & 0x3A) | ((r >> 2) & 1) | ((r << 2) & 4); break;
    }
    return *this;
}
const char* OrthoRef3::oristring ="\x00\x01\x02\x00\x01\x00\x02\x01\x00";
const char* OrthoRef3::getOri() const{
	switch(r >> 3){
		case 1: return oristring + 5;
		case 2: return oristring + 1;
		case 3: return oristring + 4;
		case 4: return oristring + 2;
		case 5: return oristring + 6;
		default: return oristring;
	}
}
const OrthoRef3& OrthoRef3::show(FILE*f, int level) const{
    unsigned short dir = (r >> 4) & 3;
    char chout[13];
    printf("%i\n", r);
    if (dir == 3) fprintf(f,"Not a direction:\n000\n000\n000\n");
    else{
    memcpy(chout, "000\n000\n000\n", sizeof(char) *13);
    chout[(dir % 3) << 2] = (r & 1) ? 'N' : 'P';
    chout[(((dir+ ((r & 8) ? 2 : 1)) % 3)<<2 ) |1] = (r & 2) ? 'N' : 'P';
    chout[(((dir+ ((r & 8) ? 1 : 2)) % 3)<<2 ) |2] = (r & 4) ? 'N' : 'P';
    fprintf(f,"direction:\n%s",chout);
    }
return *this;}
ratio::ratio(double val){
	vector<int> fser;
	fser.push_back((int)val);
	den = 1;
	int n;
	double r = (val - num);
	while(r != 0.0f){
		r = 1.0f / r;
		n = (int) r;
		if (((0xFFFFFFFF / n) < den)||(n > 255)) break;
		fser.push_back(n);
		den *=n;
		r -= n;
	}
	n = fser.size()-1;
	num = fser[n];
	den = 1;
	int t;
	for(n--;n>=0;n--){t = num; num = fser[n]*num + den; den = t;}
}

double timestep = 1.0f;
// LinkRegister LR;

SETCMP_enum gagahaga(const KeyElem<HyperPosition<uint32_t,3,0>, uint32_t> &mi,const KeyElem<HyperPosition<uint32_t,3,0>, uint32_t> &ma ){return(SETCMP_EQUAL);}

template < > mantissa<char,true>::operator double(){ return(((double)data) * (1.0f / 127.0f)  );}
template < > mantissa<unsigned char,true>::operator double(){ return(((double)data) * (1.0f / 255.0f) );}
template < > mantissa<short,true>::operator double(){ return(((double)data) * (1.0f / 32767.0f) );}
template < > mantissa<unsigned short,true>::operator double(){return(((double)data) * (1.0f / 65535.0f) );}
template < > mantissa<int,true>::operator double(){return(((double)data) * (1.0f / 2147483647.0f) );}
template < > mantissa<uint32_t,true>::operator double(){return(((double)data) * (1.0f / 4294967295.0f) );}

template < > mantissa<char,false>::operator double(){ return(((double)data) * pow(0.5f, 7.0f));}
template < > mantissa<unsigned char,false>::operator double(){ return(((double)data) * pow(0.5f, 8.0f) );}
template < > mantissa<short,false>::operator double(){ return(((double)data) * pow(0.5f, 15.0f) );}
template < > mantissa<unsigned short,false>::operator double(){return(((double)data) * pow(0.5f, 16.0f) );}
template < > mantissa<int,false>::operator double(){return(((double)data) * pow(0.5f, 31.0f) );}
template < > mantissa<uint32_t,false>::operator double(){return(((double)data) * pow(0.5f, 32.0f) );}

template < > mantissa<char,true>& mantissa<char,true>::operator=(double const & v){ data = (char)((v * 127.0f) + 0.5f);return(*this);}
template < > mantissa<unsigned char,true>& mantissa<unsigned char,true>::operator=(double const & v){data = (unsigned char)((v * 255.0f) + 0.5f);return(*this);}
template < > mantissa<short,true>& mantissa<short,true>::operator=(double const & v){data = (short)((v * 32767.0f) + 0.5f);return(*this);}
template < > mantissa<unsigned short,true>& mantissa<unsigned short,true>::operator=(double const & v){data = (unsigned short)((v * 65535.0f) + 0.5f);return(*this);}
template < > mantissa<int,true>& mantissa<int,true>::operator=(double const & v){data = (int)((v * 2147483647.0f) + 0.5f);return(*this);}
template < > mantissa<uint32_t,true>& mantissa<uint32_t,true>::operator=(double const & v){data = (uint32_t)((v * 4294967295.0f) + 0.5f);return(*this);}

template < > mantissa<char,false>& mantissa<char,false>::operator=(double const & v){data = (char)((v * 128.0f) + 0.5f);return(*this);}
template < > mantissa<unsigned char,false>& mantissa<unsigned char,false>::operator=(double const & v){data = (unsigned char)((v * 256.0f) + 0.5f);return(*this);}
template < > mantissa<short,false>& mantissa<short,false>::operator=(double const & v){data = (short)((v * 32768.0f) + 0.5f);return(*this);}
template < > mantissa<unsigned short,false>& mantissa<unsigned short,false>::operator=(double const & v){data = (unsigned short)((v * 65536.0f) + 0.5f);return(*this);}
template < > mantissa<int,false>& mantissa<int,false>::operator=(double const & v){data = (int)((v * 2147483648.0f) + 0.5f);return(*this);}
template < > mantissa<uint32_t,false>& mantissa<uint32_t,false>::operator=(double const & v){data = (uint32_t)((v * 4294967296.0f) + 0.5f);return(*this);}

Tuple<double, 3> HCYfromRGB(Tuple<double, 3> rgb, bool integer_hue){ Tuple<double, 3> fout;
    double m;
    fout[1] = 1.0f;
   // fout[2] = rgb[0] * 0.21875f + rgb[1] * 0.6875f + rgb[2] * 0.09375f;
    fout[2] = rgb[0] * 0.375f + rgb[1] * 0.5f + rgb[2] * 0.125f;

    if (rgb[0] >= rgb[1]){
        if (rgb[1] >= rgb[2]){ // RGB 0
            m = rgb[0] - rgb[2];
            if (m == 0.0f) { // grayscale
                fout[0] = 0;
                fout[1] = 0;
            }else{
                fout[0] = (rgb[1] - rgb[2]) / m;
                if (fout[2] < 0.375f + 0.5f * fout[0]){
                    fout[1] = 1.0f - rgb[2] / fout[2];
                }else{
                    fout[1] = (rgb[0] - fout[2]) / (1.0f - fout[2]);
                }
            }
        }else if (rgb[0] > rgb[2]){ // RBG 5
            m = rgb[0] - rgb[1];
            fout[0] = (rgb[2] - rgb[1]) / m;  //(hm)
            if (fout[2] < 0.375f + 0.125f * fout[0]){
                fout[1] = 1.0f - rgb[1] / fout[2];
            }else{
                fout[1] = (rgb[0] - fout[2]) / (1.0f - fout[2]);
            }
            fout[0] = 6.0f - fout[0];
        }else{ // BRG 4
            m = rgb[2] - rgb[1];
            fout[0] = (rgb[0] - rgb[1]) / m;  //(hm)
            if (fout[2] < 0.125f + 0.375f * fout[0]){
                fout[1] = 1.0f - rgb[1] / fout[2];
            }else{
                fout[1] = (rgb[2] - fout[2]) / (1.0f - fout[2]);
            }
            fout[0] += 4.0f;
        }
    }else{
        if (rgb[0] > rgb[2]){ // GRB 1
            m = rgb[1] - rgb[2];
            fout[0] = (rgb[0] - rgb[2]) / m;  //(hm)
            if (fout[2] < 0.5f + 0.375f * fout[0]){
                fout[1] = 1.0f - rgb[2] / fout[2];
            }else{
                fout[1] = (rgb[1] - fout[2]) / (1.0f - fout[2]);
            }

            fout[0] = 2.0f - fout[0];
        }else if (rgb[1] > rgb[2]){ // GBR 2
            m = rgb[1] - rgb[0];
            fout[0] = (rgb[2] - rgb[0]) / m;  //(hm)

            if (fout[2] < 0.5f + 0.125f * fout[0]){
                fout[1] = 1.0f - rgb[0] / fout[2];
            }else{
                fout[1] = (rgb[1] - fout[2]) / (1.0f - fout[2]);
            }

            fout[0] += 2.0f;
        }else{ // BGR // 3
            m = rgb[2] - rgb[0];
            fout[0] = (rgb[1] - rgb[0]) / m;  //(hm)
            if (fout[2] < 0.125f + 0.5f * fout[0]){
                fout[1] = 1.0f - rgb[0] / fout[2];
            }else{
                fout[1] = (rgb[2] - fout[2]) / (1.0f - fout[2]);
            }
            fout[0] = 4.0f - fout[0];
        }
    }
    if (!integer_hue) fout[0] *= (M_PI / 3.0f);
    return fout;
}

 //   fout[2] = rgb[0] * 0.21875f + rgb[1] * 0.6875f + rgb[2] * 0.09375f;

Tuple<double, 3> RGBfromHCY(double h, double c, double y, bool integer_hue){ Tuple<double, 3> fout;
    if (!integer_hue) h /= (M_PI / 3.0f);
    float hm = fmod(h,1);
    int phase = ((int)floor(h)) % 6;
    if (h < 0) {hm += 1.0f; phase += 6;}
    float d;
    fout[0] = y * (1- c);
    fout[1] = fout[0];
    fout[2] = fout[0];
    Tuple<double, 3> query;

    //double ly = pow(y,);

    switch(phase){
        default:
            if (h > 7){
            fout[0] = y;
            fout[1] = y;
            fout[2] = y;
            break;
            }
        case 0:
        if (y < 0.375f + 0.5f * hm){
            d = y / (0.375f + 0.5f * hm);
            fout[0] += c* d;
            fout[1] += c* hm * d;
        }else{
            d = (1.0f-y) / (0.625f - 0.5f * hm);
            fout[0] += c;
            fout[1] += c * (1.0 - d * (1.0f- hm));
            fout[2] += c * (1.0f-d);
        }
        break;
        case 1:
        if (y < 0.875f - 0.375f * hm){
            d = y / (0.875f - 0.375f * hm);
            fout[0] += c *(1.0f-hm) * d;
            fout[1] += c *d;
            //fout[2] += 0.0f;
        }else{
            d = (1.0f-y) / (0.125f + 0.375f * hm);
            fout[0] += c*(1.0 - d * hm);
            fout[1] += c;
            fout[2] += c*(1.0f-d);
        }
/*
        if (y < 0.21875f + (0.6875f - 0.21875f) * hm){ //
            d = y / (0.21875f + (0.6875f - 0.21875f) * hm);
            fout[0] += c* d;
            fout[1] += c* hm * d;
        }else{
            d = (1.0f-y) / ( 1.0f - 0.21875f - (0.6875f - 0.21875f) * hm);
            fout[0] += c;
            fout[1] += c * (1.0 - d * (1.0f- hm));
            fout[2] += c * (1.0f-d);
        }
        break;
        case 1:
        if (y < (0.6875f + 0.21875f) - 0.21875f * hm){
            d = y / ( (0.6875f + 0.21875f) - 0.21875f * hm);
            fout[0] += c *(1.0f-hm) * d;
            fout[1] += c *d;
            //fout[2] += 0.0f;
        }else{
            d = (1.0f-y) / (1.0f - (0.6875f + 0.21875f) + 0.21875f * hm);
            fout[0] += c*(1.0 - d * hm);
            fout[1] += c;
            fout[2] += c*(1.0f-d);
        }*/
        break;
        case 2:
        if (y < 0.5f + 0.125f * hm){
            d = y / (0.5f + 0.125f * hm);
            // fout[0] += 0;
            fout[1] += c* d;
            fout[2] += c*(hm * d);
        }else{
            d = (1.0f-y) / ( 0.5f - 0.125f * hm);
            fout[0] += c * (1.0f-d);
            fout[1] += c;
            fout[2] += c * (1.0 - d * (1.0f- hm));
        }
        break;
        case 3:
        if (y < 0.625f - 0.5f * hm){
            d = y / (0.625f - 0.5f * hm);
            //fout[0] += 0.0f;
            fout[1] += c *(1.0f-hm) * d;
            fout[2] += c *d;

        }else{
            d = (1.0f-y) / (0.375f + 0.5f * hm);
            fout[0] += c*(1.0f-d);
            fout[1] += c*(1.0 - d * hm);
            fout[2] += c;
        }
        break;
        case 4:
        if (y < 0.125f + 0.375f * hm){
            d = y / (0.125f + 0.375f * hm);
            // fout[0] += 0;
            fout[2] += c* d;
            fout[0] += c*(hm * d);
        }else{
            d = (1.0f-y) / ( 0.875f - 0.375f * hm);
            fout[1] += c * (1.0f-d);
            fout[2] += c;
            fout[0] += c * (1.0 - d * (1.0f- hm));
        }
        break;
        case 5:
        if (y < 0.5f - 0.125f * hm){
            d = y / (0.5f - 0.125f * hm);
            //fout[0] += 0.0f;
            fout[2] += c *(1.0f-hm) * d;
            fout[0] += c *d;

        }else{
            d = (1.0f-y) / (0.5f + 0.125f * hm);
            fout[1] += c*(1.0f-d);
            fout[2] += c*(1.0 - d * hm);
            fout[0] += c;
        }
        break;
    }
    return fout;
}



Tuple<double, 3> RGBfromHCY_gamma(double h, double c, double y, double g, bool opt_hue){ Tuple<double, 3> fout;

    if (opt_hue) h = h - sin(h *18.8495559215) / 24;
    y = pow(y,g);

    double res = (h < 0.0) ? (1.0 - fmod(-h, 1)) * 6 : fmod(h, 1) * 6;
    int bas = (int) floor(res);
    res = res - bas;
    switch(bas){
    case 0:
        fout[0] = 1.0; fout[1] = 1.0 - c * (1-res); fout[2] = 1.0 - c;
    break; case 1:
        fout[0] = 1.0 - c * res; fout[1] = 1.0; fout[2] = 1.0 - c;
    break; case 2:
        fout[0] = 1.0 - c; fout[1] = 1.0; fout[2] = 1.0 - c* (1-res);
    break; case 3:
        fout[0] = 1.0 - c; fout[1] = 1.0 - c*res; fout[2] = 1.0;
    break; case 4:
        fout[0] = 1.0 - c * (1-res); fout[1] = 1.0 - c; fout[2] = 1.0;
    break; case 5:
        fout[0] = 1.0; fout[1] = 1.0 - c; fout[2] = 1.0-c*res;
    }
    double brscale = fout[0] * 0.299 + fout[1] * 0.587 + fout[2] * 0.114;
    double deriv,tmp, discr;
    if (brscale < y){
        if ((fout[0] < 1)&(fout[0] >0.00001)) {tmp = fout[0]; fout[0] = pow(fout[0],(1/g)) ; deriv = 0.3*(1/fout[0]-1)*tmp; }
        else deriv =0;
        if ((fout[1] < 1)&(fout[1] >0.00001)) {tmp = fout[1]; fout[1] = pow(fout[1],(1/g)) ; deriv += 0.59*(1/fout[1]-1)*tmp;}
        if ((fout[2] < 1)&(fout[2] >0.00001)) {tmp = fout[2]; fout[2] = pow(fout[2],(1/g)) ; deriv += 0.11*(1/fout[2]-1)*tmp;}

        deriv = deriv * g;
        discr = deriv * deriv - 4.0 * (brscale -y) * (1.0 -brscale - deriv);
        brscale = (-deriv + sqrt(discr)) / (2.0* (1.0 -brscale - deriv));
        fout[0] = fout[0] + (1 -fout[0]) * brscale;
        fout[1] = fout[1] + (1 -fout[1]) * brscale;
        fout[2] = fout[2] + (1 -fout[2]) * brscale;
    }else{
        brscale =  y / brscale;
        g = 1.0 / g;
        fout[0] = pow((fout[0]* brscale), g);
        fout[1] = pow((fout[1]* brscale), g);
        fout[2] = pow((fout[2]* brscale), g);
    }
return fout;}

Tuple<double, 3> RGBfromHCY_gamma2(double h, double c, double y, bool opt_hue){ Tuple<double, 3> fout;

    if (opt_hue) h = h - sin(h *18.8495559215) / 24;

    y = (y < 0.0392768) ? y / 12.92 : pow((y + 0.055)/ 1.055, 2.4);


    double res = fmod(h, 1) * 6;
    int bas = (int) floor(res);
    res = res - bas;
    switch(bas){
    case 0:
        fout[0] = 1.0; fout[1] = 1.0 - c * (1-res); fout[2] = 1.0 - c;
    break; case 1:
        fout[0] = 1.0 - c * res; fout[1] = 1.0; fout[2] = 1.0 - c;
    break; case 2:
        fout[0] = 1.0 - c; fout[1] = 1.0; fout[2] = 1.0 - c* (1-res);
    break; case 3:
        fout[0] = 1.0 - c; fout[1] = 1.0 - c*res; fout[2] = 1.0;
    break; case 4:
        fout[0] = 1.0 - c * (1-res); fout[1] = 1.0 - c; fout[2] = 1.0;
    break; case 5:
        fout[0] = 1.0; fout[1] = 1.0 - c; fout[2] = 1.0-c*res;
    }
    double brscale = fout[0] * 0.299 + fout[1] * 0.587 + fout[2] * 0.114;
    double deriv,tmp, discr;
    if (brscale < y){
        if ((fout[0] < 1)&&(fout[0] >0.00001)) {tmp = fout[0]; fout[0] = (fout[0] < 0.00304) ? 12.92 * fout[0] : (1.055 * pow(fout[0], 1.0 / 2.4) - 0.055) ; deriv = 0.3*(1/fout[0]-1)*tmp; }
        else deriv =0;
        if ((fout[1] < 1)&&(fout[1] >0.00001)) {tmp = fout[1]; fout[1] = (fout[1] < 0.00304) ? 12.92 * fout[1] : (1.055 * pow(fout[1], 1.0 / 2.4) - 0.055) ; deriv += 0.59*(1/fout[1]-1)*tmp;}
        if ((fout[2] < 1)&&(fout[2] >0.00001)) {tmp = fout[2]; fout[2] = (fout[2] < 0.00304) ? 12.92 * fout[2] : (1.055 * pow(fout[2], 1.0 / 2.4) - 0.055) ; deriv += 0.11*(1/fout[2]-1)*tmp;}
        deriv = deriv * 2.4;
        discr = deriv * deriv - 4.0 * (brscale -y) * (1.0 -brscale - deriv);
        brscale = (-deriv + sqrt(discr)) / (2.0* (1.0 -brscale - deriv));
        fout[0] = fout[0] + (1 -fout[0]) * brscale;
        fout[1] = fout[1] + (1 -fout[1]) * brscale;
        fout[2] = fout[2] + (1 -fout[2]) * brscale;
    }else{
        brscale =  y / brscale;
        fout[0] *= brscale; fout[1] *= brscale; fout[2] *= brscale;
        fout[0] = (fout[0] < 0.00304) ? 12.92 * fout[0] : (1.055 * pow(fout[0], 1.0 / 2.4) - 0.055); // pow((fout[0]* brscale), g);
        fout[1] = (fout[1] < 0.00304) ? 12.92 * fout[1] : (1.055 * pow(fout[1], 1.0 / 2.4) - 0.055);
        fout[2] = (fout[2] < 0.00304) ? 12.92 * fout[2] : (1.055 * pow(fout[2], 1.0 / 2.4) - 0.055);
    }
return fout;}

double gammaCorrHue(double hue, double gamma){
    int i = floor(hue);
    if (i & 1) return 1.0f - pow(1.0f - (hue - i), gamma) + i;
    else return pow(hue - i, gamma) + i;
}

/*
Tuple<double, 3> RGBfromHCY(double h, double c, double y, bool integer_hue){ Tuple<double, 3> fout;
    if (!integer_hue) fout[0] /= (M_PI / 3.0f);
    float hm = fmod(h,1);
    int phase = ((int)floor(h)) % 6;
    if (h < 0) {hm += 1.0f; phase += 6;}
    float d;
    fout[0] = y * (1- c);
    fout[1] = fout[0];
    fout[2] = fout[0];
    Tuple<double, 3> query;
    switch(phase){
        default:
            if (h > 7){
            fout[0] = y;
            fout[1] = y;
            fout[2] = y;
            break;
            }
        case 0:
        if (y < 0.375f + 0.5f * hm){
            d = y / (0.375f + 0.5f * hm);
            fout[0] += c* d;
            fout[1] += c* hm * d;
        }else{
            d = (1.0f-y) / (0.625f - 0.5f * hm);
            fout[0] += c;
            fout[1] += c * (1.0 - d * (1.0f- hm));
            fout[2] += c * (1.0f-d);
        }
        break;
        case 1:
        if (y < 0.875f - 0.375f * hm){
            d = y / (0.875f - 0.375f * hm);
            fout[0] += c *(1.0f-hm) * d;
            fout[1] += c *d;
            //fout[2] += 0.0f;
        }else{
            d = (1.0f-y) / (0.125f + 0.375f * hm);
            fout[0] += c*(1.0 - d * hm);
            fout[1] += c;
            fout[2] += c*(1.0f-d);
        }
        break;
        case 2:
        if (y < 0.5f + 0.125f * hm){
            d = y / (0.5f + 0.125f * hm);
            // fout[0] += 0;
            fout[1] += c* d;
            fout[2] += c*(hm * d);
        }else{
            d = (1.0f-y) / ( 0.5f - 0.125f * hm);
            fout[0] += c * (1.0f-d);
            fout[1] += c;
            fout[2] += c * (1.0 - d * (1.0f- hm));
        }
        break;
        case 3:
        if (y < 0.625f - 0.5f * hm){
            d = y / (0.625f - 0.5f * hm);
            //fout[0] += 0.0f;
            fout[1] += c *(1.0f-hm) * d;
            fout[2] += c *d;

        }else{
            d = (1.0f-y) / (0.375f + 0.5f * hm);
            fout[0] += c*(1.0f-d);
            fout[1] += c*(1.0 - d * hm);
            fout[2] += c;
        }
        break;
        case 4:
        if (y < 0.125f + 0.375f * hm){
            d = y / (0.125f + 0.375f * hm);
            // fout[0] += 0;
            fout[2] += c* d;
            fout[0] += c*(hm * d);
        }else{
            d = (1.0f-y) / ( 0.875f - 0.375f * hm);
            fout[1] += c * (1.0f-d);
            fout[2] += c;
            fout[0] += c * (1.0 - d * (1.0f- hm));
        }
        break;
        case 5:
        if (y < 0.5f - 0.125f * hm){
            d = y / (0.5f - 0.125f * hm);
            //fout[0] += 0.0f;
            fout[2] += c *(1.0f-hm) * d;
            fout[0] += c *d;

        }else{
            d = (1.0f-y) / (0.5f + 0.125f * hm);
            fout[1] += c*(1.0f-d);
            fout[2] += c*(1.0 - d * hm);
            fout[0] += c;
        }
        break;
    }
    return fout;
}*/

void RGB2Dvec(unsigned char *fout, double *f_in){

    double norm = sqrt(f_in[0]*f_in[0] + f_in[1] *f_in[1]);
    if (norm == 0.0f){
    fout[0] = 0;
    fout[1] = 0;
    fout[2] = 0;
    }else if (norm < 1.0f){
    double unitgamma = pow(-f_in[0]+norm, 2.2f) + 1.9667f * pow(0.1499999f*f_in[1]+ 0.20f * f_in[0]+ norm * 0.5f , 2.2f) + 0.36f * pow(-0.5999999f*f_in[0]+ 0.8f * f_in[1] + norm, 2.2f);
    // real gamma = pow(F,2.2)
    // output gamma = 0.3 8 (exp^- norm) -1.0f

    unitgamma = 295.0f * pow((1.0f - exp(-0.4375f *norm)) / unitgamma, 1.0f / 2.2f);

    fout[0] = (unsigned char) (unitgamma * (-f_in[0]+norm));
    fout[1] = (unsigned char) (unitgamma * (0.1499999f*f_in[1]+ 0.20f * f_in[0]+ norm * 0.5f ));
    fout[2] = (unsigned char) (unitgamma * ( f_in[1]+ norm));

    }else{
    double unitgamma = pow(-f_in[0]+norm, 2.2f) + 1.9667f * pow(0.1499999f*f_in[1]+ 0.20f* f_in[0]+ 0.5f *norm * (2.0f - exp(1.0f-norm)) , 2.2f) + 0.36f * pow(-0.5999999f*f_in[0]+ 0.8f * f_in[1] + norm, 2.2f);
    // real gamma = pow(F,2.2)
    // output gamma = 0.3 8 (exp^- norm) -1.0f

    unitgamma = 295.0f * pow((1.0f - exp(-0.4375f *norm)) / unitgamma, 1.0f / 2.2f);

    fout[0] = (unsigned char) (unitgamma * (-f_in[0]+norm));
    fout[1] = (unsigned char) (unitgamma * (0.1499999f*f_in[1]+ 0.20f * f_in[0]+ 0.5 * norm *(2.0f - exp(1.0f-norm))));
    fout[2] = (unsigned char) (unitgamma * (-0.5999999f*f_in[0]+ 0.8f * f_in[1]+ norm));
    }
  // 1.0 red, 0.75 green 1.0 blue
  // objective : 0.3 after 2.2 gamma
  // red = F * (-f_in[0] + norm) * 0.5f
  // green = F * (f_in[0] + norm) * 0.325f
  // blue = F * (f_in[1] + norm) * 0.5f

  // end gamma = F^2.2  * (-f_in[0] + norm) * 0.5f

}


Semaphore::ReadAccess::operator bool(){
    int preval = target.sema.fetch_add(1);
    if (preval < 256) {state =1; return true;}
    else{target.sema.fetch_add(-1); state =0; return false;}
}
Semaphore::WriteAccess::operator bool(){
    int preval = target.sema.fetch_add(256);
    if (preval == 256) {state =256; return true;}
    else{target.sema.fetch_add(-256); state =0; return false;}
}


Semaphore::WaitForReadAccess::operator bool(){
    int preval;
    while(true){
        preval = target.sema.fetch_add(1);
        if (preval < 256) break;
        target.sema.fetch_add(-1);
    }
return true;}
Semaphore::WaitForWriteAccess::operator bool(){
    int preval;
    while(true){
        preval = target.sema.fetch_add(256);
        if (preval == 0) break;
        target.sema.fetch_add(-256);
    }
return true;}

MutexAccess::ReadAccess::operator bool(){
    if (target.mtx.try_lock()){
        target.sema.fetch_add(1);
    }else{
        std::unique_lock<std::mutex> lck(target.mtx);
        target.cv.wait(lck);
    }
return true;}

MutexAccess::WriteAccess::WriteAccess(MutexAccess& _target) :lock(_target.mtx, std::defer_lock){
    int nbread = _target.sema.fetch_add(1);
    if (nbread == 0){
        nbread--;
        lock.lock();
    }else{
        // TODO
    }

}

//trajecion/ umap
//nuclei
//damethods
//velocity

ScriptCompiler::ScriptCompiler():errbuffer(NULL), nbdeclared(0){
	types.createEntryWithAlias("byte", 0x0100);
	types.createEntryWithAlias("short", 0x0101);
	types.createEntryWithAlias("int", 0x0102);
	types.createEntryWithAlias("long", 0x0103);
	types.createEntryWithAlias("void", 0x0104);
	types.createEntryWithAlias("bool", 0x0105);
	types.createEntryWithAlias("float", 0x0106);
	types.createEntryWithAlias("double", 0x0107);
}
void ScriptCompiler::initFnccodes(){
	Tuple<uint16_t> fnccodes_ins;
	fnccodes_ins.setSize(2);
	fnccodes_ins[0] = 0x0107;	fnccodes_ins[1] = 0x0107;
	fnccodes.createEntryWithAlias(string("sqrt"), 0xFFFE) = fnccodes_ins;
	fnccodes.createEntryWithAlias(string("sin"), 0xFFFD) = fnccodes_ins;
	fnccodes.createEntryWithAlias(string("cos"), 0xFFFC) = fnccodes_ins;
}

uint16_t ScriptCompiler::addEnumType(string type){

return 0;}
void addEnumType(uint16_t enumCode, const char* type, uint32_t value){

}
void ScriptCompiler::declareFunction(string name, Tuple<const char*> arg_types ){
	Tuple<uint16_t> xf_arg_types; xf_arg_types.setSize( arg_types.getSize());
	for(unsigned int i=0;i< arg_types.getSize();i++) xf_arg_types[i] = this->getTypeCode(arg_types[i]);
return declareFunction(name,xf_arg_types);}
void ScriptCompiler::declareFunction(string name, Tuple<uint16_t> arg_types ){
	nbdeclared++;
	fnccodes.createEntryWithAlias(name, nbdeclared) =arg_types;
return;}
void ScriptCompiler::defineFunction(string name, std::function< char* (Ambiguous*&) > dafnc ){
	uint32_t func = fnccodes.getAlias(name);
	if (func == 0) {printf("Trying to define function %s which was not declared!\n", name.c_str()); exit(1);}
	if (func >= externfnc.getSize()) externfnc.upSize(func+1);
	externfnc[func] = dafnc;
return;}


/*
void ScriptCompiler::addFunction(string name, Tuple<uint16_t> arg_types , std::function< void (Ambiguous*&) > dafnc ){

}

void ScriptCompiler::addExternFunction(string name, Tuple<const char*> arg_types , uint32_t fnccode ){
	Tuple<uint16_t> xf_arg_types; xf_arg_types.setSize( arg_types.getSize());
	for(int i=0;i< arg_types.getSize();i++) xf_arg_types[i] = this->getTypeCode(arg_types[i]);
return addExternFunction(name,xf_arg_types,fnccode);}
void ScriptCompiler::addExternFunction(string name, Tuple<uint16_t> arg_types , uint32_t fnccode ){
}*/

int ScriptCompiler::mkCommand_matchtypes_routine(StructuredHash<uint16_t, uint32_t>::TreeIterator<0> &ite, uint32_t dacode, uint32_t trmax){
	auto nite = ite + 1;
	//uint32_t ntype = (*nite);
	nite++ ; *nite = dacode;
	if ((*ite) == (*nite)) {
		if (ite() > trmax) *ite = 0;
		ite++;
		if (ite() > trmax) *ite = 0;
	}else switch(*ite){
		case 0x0104: // logical, these operators are illegal, but silly conversions? maybe?
			printf("warning, trying operator on logical value\n");

		case 0x0100: case 0x0101: case 0x0102: case 0x0103:
		if (ite() > trmax) *ite = 0;
		ite++;
		switch(*nite){
			case 0x0106:
			if (ite() > trmax) *ite = LFHSCRIPT_64BITS_TO_F32BITS;
			else return 1;
		break; case 0x0107:
			if (ite() > trmax) *ite = LFHSCRIPT_64BITS_TO_F64BITS;
			else return 1;
		break; default:
			if (ite() > trmax) *ite = 0; // all good
		break;}
	break; case 0x0106:
		if (ite() > trmax) *ite = 0;
		ite++;
		switch(*nite){
			case 0x0107:
			if (ite() > trmax) *ite = LFHSCRIPT_F64BITS_TO_F32BITS;
		break; default:
			if (ite() > trmax) *ite = LFHSCRIPT_64BITS_TO_F32BITS;
			else return 1;
		break;}
	break; case 0x0107:
		if (ite() > trmax) *ite = 0;
		ite++;
		switch(*nite){
			case 0x0106:
			if (ite() > trmax) *ite = LFHSCRIPT_F32BITS_TO_F64BITS;
			else return 1;
		break; default:
			if (ite() > trmax) *ite = LFHSCRIPT_64BITS_TO_F64BITS;
			else return 1;
		}
	break;
	}
return 0;}

bool ScriptCompiler::printError(FILE* f){
	if (errbuffer == NULL) return false;
	fprintf(f, "%s", errbuffer);
	delete[](errbuffer); errbuffer = NULL;
return true;}

void ScriptCompiler::errorSet(const char * fnc){
	if (errbuffer == NULL) {
		int i = strlen(fnc)+1;
		errbuffer = new char[i];
		memcpy(errbuffer, fnc, i);
	}
}
void ScriptCompiler::errorSet(const char * fnc, char dachar){
	if (errbuffer == NULL) {
		int i = strlen(fnc)+8;
		errbuffer = new char[i];
		sprintf(errbuffer, fnc, dachar);
	}
}
void ScriptCompiler::errorSet(const char * fnc, const char *dastr){
	if (errbuffer == NULL) {
		int i = strlen(fnc) + strlen(dastr)+8;
		errbuffer = new char[i];
		sprintf(errbuffer, fnc, dastr);
	}
}


void ScriptCompiler::errorNoArgumentMatch_routine(string fnc, const char* precise){
	uint32_t i;
	if (errbuffer == NULL) errbuffer = new char[2048];
	if (auto ite = fnccodes.mkIterator(fnc.c_str())) {
		uint32_t hehe = sprintf(errbuffer, "Function '%s' could not find a match %s: Valid candidates are:\n", fnc.c_str(), precise);
		do{
			hehe += sprintf(errbuffer + hehe, "\t%s %s", types.deref_key(U32_Alias((*ite)[0])).c_str(), fnc.c_str());
			for(i=1;i< ite->getSize();i++){
				hehe += sprintf(errbuffer + hehe, "%c%s%c", (i == 1) ? '(' : ' ',  types.deref_key(U32_Alias((*ite)[i])).c_str(), i== ite->getSize() -1 ? ')': ',');
			}
			hehe += sprintf(errbuffer + hehe, (i == 1) ? "()\n" : "\n");
		}while(ite++);
	}else {
		sprintf(errbuffer, "'%s' is not a known function\n", fnc.c_str());
	}
}

Script ScriptCompiler::mkCommand(const char* command){ Script fout(*this);
	myHashmap<uint32_t, Tuple<uint32_t, 2> > parser;
	myHashmap<uint32_t,  Ambiguous > constants;

	//printf("Trying to compile: %s\n", command);

	myHashmap<uint32_t, KeyElem<uint32_t, uint16_t> > symbolbuf;
	StructuredHash<uint16_t, uint32_t> regular;  // id, symbol

	uint8_t priority[4096];
	// Expr -> (Expr)
	// Expr -> Expr Oper Expr
	// Expr -> Func(ExprCommaList)
	// ExprCommaList -> ""
	// ExprCommaList -> Expr,ExprCommaList
	// Expr -> Expr ? Expr : Expr
	// Expr -> constant
	// constant -> - constant
	// Oper -> + - * /
	// Func -> sin cos ...

	uint32_t i,j;
	Vector<uint32_t> trace;
	Vector<uint16_t> strace;
	const char* cur = command;
	const char* scur;
	string astring;

	// variable types: logical int float enum
	// 8bits 16bits 32bits 64bits logical na.. float double

	// step 1: parse literals and lookup function names
	while(*cur != '\0'){
		switch(*cur){
		case '\n':
		case '\r':
		case '\t':
		case ' ':
			cur++;
		break; case '!': if (cur[1] != '=') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '"': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '#': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '$': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '%': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '&': if (cur[1] != '&') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '\'': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '(': strace.push_back((uint16_t)*(cur++)); // symbolbuf[symbufcur] = KeyElem<uint32_t, uint16_t>(symbufcur==0?1:2,*(cur++)); symbufcur++;
		break; case ')': strace.push_back((uint16_t)*(cur++)); // symbolbuf[symbufcur] = KeyElem<uint32_t, uint16_t>(symbufcur==0?1:2, *(cur++)); symbufcur++;
		break; case '*': strace.push_back((uint16_t)*(cur++));
		break; case '+': strace.push_back((uint16_t)*(cur++));
		break; case ',': strace.push_back((uint16_t)*(cur++)); // symbolbuf[symbufcur] = KeyElem<uint32_t, uint16_t>(symbufcur==0?1:2, *(cur++)); symbufcur++;
		break; case '-': strace.push_back((uint16_t)*(cur++));
		break; case '.': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '/': strace.push_back((uint16_t)*(cur++));
		break; case ':': strace.push_back((uint16_t)*(cur++));
		break; case ';': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '<': if (cur[1] != '=') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '=': if (cur[1] != '=') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '>': if (cur[1] != '=') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '?': strace.push_back((uint16_t)*(cur++));
		break; case '@': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '[': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '\\': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case ']': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '^': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '`': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '{': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; case '|': if (cur[1] != '|') strace.push_back((uint16_t)*(cur++)); else {strace.push_back(0x80 | (uint16_t)*(cur++));cur++;}
		break; case '}': errorSet("%c, is not recognized character!\n", *cur); return(fout);
	//	break; case '_': errorSet("%c, is not recognized character!\n", *cur); return(fout);
		break; default:
			if (((*cur) & 64) != 0){// letter or '_'
				scur = cur +1;
				strace.push_back((uint16_t)'F');
				while (((*scur) >= '0')&&((*scur) <= 'z')&&(((*scur) >= 'a')||((*scur) <= '9')||(((*scur) >= 'A')&&(((*scur) <= 'Z')||(*scur) == '_')))) scur++;
				//constants.last().d.ptr = new char[(int)(scur - cur)];
				//memcpy(constants.last().d.ptr, cur, ((int)(scur - cur))-1);
				//((char*)constants.last().d.ptr)[((int)(scur - cur))-1] = '\0';
				astring = string(cur, ((int)(scur - cur)));
				if (auto ite = fnccodes.mkIterator(astring)) {
					//constants.last().k = LFH_ANYTHING_STRING;
					constants[strace.getSize()].int32_val = ite.getAlias();
//					printf("ok, function %s should C[%i] = %i\n", astring.c_str(), strace.getSize(), ite.getAlias());
				}else {errorSet("'%s' is not a known command!",astring.c_str()); return(fout);}

				cur = scur;
			}else{// number
				if (((*cur) == '0')&&(cur[1] == 'x')){
					// Hexadecimal

					trace.push_back((uint32_t)(cur - command));
					scur = cur+2;
					while( ((*scur >= '0')&&(*scur <= '9'))||((*scur >= 'A')&&(*scur <= 'F'))  ) scur++;
					bool tmp;
					bool isunsigned = false;
					int width = scur - cur -2;
					while( ( *scur == 'L')||( *scur == 'l')||(tmp = ((*scur == 'u')||( *scur == 'U')))) {scur++; isunsigned |= tmp;}
					int hehe;
					constants[strace.getSize()+1].int64_val = 0;
					if  (width <= 4){
						sscanf(cur, "%i", &hehe);
						if  (width <= 2){
							strace.push_back(0x0100);
							constants[strace.getSize()].int8_val = hehe;
							//constants.last().k = (isunsigned) ? LFH_ANYTHING_UNSIGNED_CHAR : LFH_ANYTHING_CHAR;
						}else{
							strace.push_back(0x0101);
							constants[strace.getSize()].int16_val = hehe;
							//constants.last().k = (isunsigned) ? LFH_ANYTHING_UNSIGNED_SHORT : LFH_ANYTHING_SHORT;
						}
					}else if  (width <= 8){
						sscanf(cur, "%i", &hehe);
						strace.push_back(0x0102);
						constants[strace.getSize()].int32_val = hehe;
						//constants.last().k = (isunsigned) ? LFH_ANYTHING_UNSIGNED_INT : LFH_ANYTHING_INT;
					}else{
						sscanf(cur, "%i", &hehe);
						strace.push_back(0x0103);
						constants[strace.getSize()].int32_val = hehe;
						//constants.last().k = (isunsigned) ? LFH_ANYTHING_UNSIGNED_LONG : LFH_ANYTHING_LONG;
					}
					cur = scur;
				}else{ // int or float

					scur = cur +1;
					while( ((*scur >= '0')&&(*scur <= '9'))) scur++;
					switch(*scur){
						case '.':
						scur++;
						while( ((*scur >= '0')&&(*scur <= '9'))) scur++;
						sscanf(cur, "%f", &constants[strace.getSize()+1].float_val);
						if ((*scur) != 'f') {strace.push_back(0x0107); constants[strace.getSize()].double_val = (double)constants[strace.getSize()].float_val; }
						else {scur++; strace.push_back(0x0106);}
					break; case 'f': scur++; strace.push_back(0x0106); sscanf(cur, "%f", &constants[strace.getSize()].float_val);
					break; default: strace.push_back(0x0102); sscanf(cur, "%i", &constants[strace.getSize()].int32_val);
					}
					cur = scur;
				}

			}
		}
	}
	regular.insertElem(0, (uint32_t)'^');
	for(uint32_t state = strace.getSize(); state>0; state--){
		regular.insertElemAfter(state, strace[state-1] , 0);
	}
	//regular.show();

	// check () [] "" '' {} matching



	uint32_t symbufcur = 1;
	Vector< HeapTree< KeyElem<uint32_t, Tuple<uint32_t, 3u> >  > > ordered_resols;
	myHashmap<uint32_t, Tuple<uint32_t, 4u> > rule_resols; // (extra, first, second, sum)
	KeyElem<uint32_t, Tuple<uint32_t, 3u> > ordered_resols_elem;


	// init Heaptree or potential resolutions

	/*
	if (auto ite = regular.mkIterator(0)) do{
		if ((*ite) < 128) printf("%i:%c\t", ite(), (char)(*ite));
		else printf("%i:%X\t", ite(), (*ite));
	}while(ite.next()); printf("\n");*/
	// scan current array for rule matches

	Vector< char > parenthesis;
	Vector< uint32_t > partial;
	Vector< KeyElem<uint32_t, uint32_t> > complete;
	symbufcur = 0;
	uint32_t par_depth = 0; // counts matching depth for () {} [] pairs



	if (auto ite = regular.mkIterator(0)) do{
		priority[symbufcur] = 0;
		for(uint32_t pcur =0; pcur < partial.getSize();pcur++){
			ordered_resols_elem.k = 0;
			switch( partial[pcur] | (((*ite) <= 0x100) ? (*ite) : 0x100)){
				case 0x00100000 | '!': partial[pcur] = 0x02210000;   // unnary !
			break; case 0x00100000 | '-': partial[pcur] = 0x022D0000;   // unnary -
			break; case 0x00100000 | '*': partial[pcur] = 0x022A0000;   // unnary *
			break; case 0x00100000 | '&': partial[pcur] = 0x02260000;   // unnary &
			break; case 0x00030100: partial.pop_swap(pcur--); // -i detected
			break; case 0x00050100: partial.pop_swap(pcur--); // !i detected
			break; case 0x00040000 | '(': partial[pcur] = 0x01040000; // f( detected
			break; case 0x01040100: partial[pcur] = 0x02040000; // f(i detected
			break; case 0x01040000 | 'c': partial[pcur] = 0x03040000; // f(c detected
			break; case 0x00020100: partial[pcur] = 0x01020000; // (i
			break; case 0x00060000 | ',': partial[pcur] = 0x01C10100; // c,

			break; case 0x00010000 | '?': partial[pcur] = 0x003F0000; // e?
			break; case 0x00010000 | '+': partial[pcur] = 0x002B0000;
			break; case 0x00010000 | '-': partial[pcur] = 0x002D0000;
			break; case 0x00010000 | '*': partial[pcur] = 0x002A0000;
			break; case 0x00010000 | '/': partial[pcur] = 0x002F0000;
			break; case 0x00010000 | '%': partial[pcur] = 0x00250000;

			break; case 0x00010000 | '<': partial[pcur] = 0x003C0000;
			break; case 0x00010000 | '=': partial[pcur] = 0x003D0000;
			break; case 0x00010000 | '>': partial[pcur] = 0x003E0000;
			break; case 0x00010080 | '<': partial[pcur] = 0x013C0000;
			break; case 0x00010080 | '=': partial[pcur] = 0x013D0000;
			break; case 0x00010080 | '>': partial[pcur] = 0x013E0000;
			break; case 0x00010000 | '&': partial[pcur] = 0x00260000;
			break; case 0x00010080 | '&': partial[pcur] = 0x01260000;
			break; case 0x00010000 | '|': partial[pcur] = 0x007C0000;
			break; case 0x00010080 | '|': partial[pcur] = 0x017C0000;
			break; case 0x00010080 | '!': partial[pcur] = 0x01210000;
			break; case 0x00010000 | ',': partial[pcur] = 0x00C10100;


			//break; case 0x00020000 | '+': partial[pcur] = 0x01020000;
			//break; case 0x00020000 | '-': partial[pcur] = 0x01020000;
			//break; case 0x00020000 | '*': partial[pcur] = 0x02020000;
			//break; case 0x00020000 | '/': partial[pcur] = 0x02020000;

			break; case 0x00050000 | 'e': partial[pcur] = 0x01050000;
			break; case 0x003F0100: partial[pcur] = 0x013F0000; // e?e
			break; case 0x013F0000 | ':': partial[pcur] = 0x023F0000; // e?e:
			break; case 0x023F0100: partial[pcur] = ordered_resols_elem.k = 0x505; // e?e:e


			break; case 0x01040000 | ')': ordered_resols_elem.k = 0x3; // f() detected
			break; case 0x02040000 | ')': ordered_resols_elem.k = 0x4; // f(i) detected
			break; case 0x03040000 | ')': ordered_resols_elem.k = 0x14; // f(c) detected
			break; case 0x01020000 | ')': ordered_resols_elem.k = 0x103; // (e) detected

			break; case 0x002A0100: ordered_resols_elem.k = 0x113;// e * e detected
			break; case 0x002F0100: ordered_resols_elem.k = 0x123;// e / e detected
			break; case 0x002B0100: ordered_resols_elem.k = 0x133;// e + e detected
			break; case 0x002D0100: ordered_resols_elem.k = 0x143;// e - e detected
			break; case 0x00250100: ordered_resols_elem.k = 0x153;// e % e detected

			break; case 0x02210100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary ! e detected
			break; case 0x022D0100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;  // unnary - e detected
			break; case 0x022A0100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary * e detected
			break; case 0x02260100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary & e detected

			break; case 0x003C0100: ordered_resols_elem.k = 0x303; // e < e detected
			break; case 0x013C0100: ordered_resols_elem.k = 0x323; // e <= e detected
			break; case 0x003D0100: ordered_resols_elem.k = 0x733; // e = e detected
			break; case 0x013D0100: ordered_resols_elem.k = 0x343; // e == e detected
			break; case 0x01210100: ordered_resols_elem.k = 0x353; // e != e detected
			break; case 0x003E0100: ordered_resols_elem.k = 0x313; // e > e detected
			break; case 0x013E0100: ordered_resols_elem.k = 0x333; // e >= e detected
			break; case 0x00260100: ordered_resols_elem.k = 0x403; // e & e detected
			break; case 0x01260100: ordered_resols_elem.k = 0x433; // e && e detected
			break; case 0x007C0100: ordered_resols_elem.k = 0x423; // e | e detected
			break; case 0x017C0100: ordered_resols_elem.k = 0x443; // e || e detected

			break; case 0x00C10100: ordered_resols_elem.k = 0x903; // e,e detected
			break; case 0x01C10100: ordered_resols_elem.k = 0x900; // c,e detected
			//break; case 0x02020000 | 'f': ordered_resols_elem.k = 0x53; // e * e detected
			break; default: partial.pop_swap(pcur--);
			}
			if (ordered_resols_elem.k != 0){
				partial.pop_swap(pcur--);
				if ((j = (ordered_resols_elem.k & 7)) == 0){ // ah... special case
					j = (ordered_resols_elem.k == 0x900) ? 3 : 4;
				}

				ordered_resols_elem.d[0] = symbufcur + 1 - j;
				ordered_resols_elem.d[1] = symbufcur;
				ordered_resols_elem.d[2] = symbufcur * j - ((j * (j -1)) >> 1);
				while (parenthesis.getSize() >= ordered_resols.getSize()) ordered_resols.push_back();
				ordered_resols[parenthesis.getSize()].insert(ordered_resols_elem);

			}
		}
		//printf("S[%i] = %c: %i partial still satisfied\n", symbufcur, (char) (*ite), partial.getSize());
		// add rules resolutions
		if ((*ite) >= 0x100) partial.push_back(0x10000); // implicit conversion to expression
		else switch(*ite){
			case '(': partial.push_back(0x20000);  partial.push_back(0x100000); parenthesis.push_back('(');
		break; case 'F': partial.push_back(0x40000);
		break; case 'c': partial.push_back(0x60000); partial.push_back(0x100000);
		break; case ')': if (parenthesis.getSize() == 0) {this->errorSet("Unexpected parenthesis closed ')'!\n"); return(fout);}
			if (parenthesis.last() != '(') {this->errorSet("Unexpected brace matching of '%c' and ')'!\n", (char) parenthesis.last()); return(fout);}
			parenthesis.pop_back();
		break; default:  partial.push_back(0x100000); // "non-value" , for unnary operators
		}
		//printf("at pos %i: %i partial, %i complete\n", symbufcur, partial.getSize(), ordered_resols.getSize());
		symbufcur++;
	}while(ite.next());

	if (parenthesis.getSize() != 0) {this->errorSet("Unexpected end, brace '%c' unmatched!\n", (char) parenthesis.last()); return(fout);}

	partial.toMemfree();

	//printf("time for tree making!\n");
	// build parse tree
	U32_Alias curalias;
	uint32_t partsum;
	int partsumbuf[16];

	Vector< KeyElem<uint32_t, Tuple<uint32_t, 3u> > > ordered_resols_buffer;
	KeyElem<uint32_t, KeyElem<uint32_t, Tuple<uint32_t, 3u> > > ordered_resols_best;

	//printf("nb par levels %i\n", ordered_resols.getSize());

	par_depth = ordered_resols.getSize() -1;
	while( ordered_resols.getSize() != 0){

		if (ordered_resols[par_depth].getSize() == 0) {ordered_resols.pop_back(); par_depth--; continue;}


		// nbsymbols is first &7, and
		// printf("%i\n", ordered_resols.top().d[0]);
		if (auto ite = regular.mkIterator(ordered_resols[par_depth].top().d[0])) {
			// printf("looking at pattern %i:", (ordered_resols.top().k & 7));
			if (ite.hasParent()) partsum = 0;
			else{
				partsum = ite();
				// printf( "[-%i-%i]", partsum, ite());

				if ((j = (ordered_resols[par_depth].top().k & 7)) == 0){ // ah... special case
					j = (ordered_resols[par_depth].top().k == 0x900) ? 3 : 4;
				}

				for(i = j-1  ; i-- >1  ;){
					if (!ite.next()) {partsum = 0; break;} // end reached... automatic fail
					// printf("%c", *ite);
					partsum += ite(); // printf( "[%i-%i]", partsum,ite());
				}
				if (i == 0){
					if (!ite.next()) partsum =0;
					else {partsum += ite(); // printf( "[-%i-%i]", partsum, ite());
						if (ite() != ordered_resols[par_depth].top().d[1]) partsum = 0; // printf("end missmatch\n");}
					}
				}
			}
		}else{ ordered_resols[par_depth].pop(); continue;}
		//printf(" partsum %i == %i ???\n", partsum, ordered_resols.top().d[2]);
		if (partsum != ordered_resols[par_depth].top().d[2]){ ordered_resols[par_depth].pop(); continue;}


		// Found 1 valid, need to find all other with equal priority!!!

		ordered_resols_best.d = ordered_resols.last().pop();
		while(!ordered_resols[par_depth].isEmpty() && (((ordered_resols[par_depth].top().k ^ ordered_resols_best.d.k) & 0xFF00) == 0)){


			if (auto ite = regular.mkIterator(ordered_resols[par_depth].top().d[0])) {
		// printf("looking at pattern %i:", (ordered_resols.top().k & 7));
			if (ite.hasParent()) partsum = 0;
			else{
				partsum = ite();
				// printf( "[-%i-%i]", partsum, ite());

				if ((j = (ordered_resols[par_depth].top().k & 7)) == 0){ // ah... special case
					j = (ordered_resols[par_depth].top().k == 0x900) ? 3 : 4;
				}

				for(i = j-1  ; i-- >1  ;){
					if (!ite.next()) {partsum = 0; break;} // end reached... automatic fail
					// printf("%c", *ite);
					partsum += ite(); // printf( "[%i-%i]", partsum,ite());
				}
				if (i == 0){
					if (!ite.next()) partsum =0;
					else {partsum += ite(); // printf( "[-%i-%i]", partsum, ite());
						if (ite() != ordered_resols[par_depth].top().d[1]) partsum = 0; // printf("end missmatch\n");}
					}
				}
			}
			}else{ordered_resols.last().pop(); continue;} // node is gone... really?

			if (partsum != ordered_resols[par_depth].top().d[2]) {ordered_resols.last().pop(); continue;}
			//printf("Found Valid equivalent!\n");
			// TODO
			if (ordered_resols_buffer.getSize() == 0) ordered_resols_best.k = regular.mkIterator(ordered_resols_best.d.d[0]).mkLastChild()();
			i = regular.mkIterator(ordered_resols[par_depth].top().d[0]).mkLastChild()();
			if (i > ordered_resols_best.k) ordered_resols_buffer.push_back(ordered_resols[par_depth].pop());
			else{
				ordered_resols_buffer.push_back( ordered_resols_best.d);
				ordered_resols_best.k = i;
				ordered_resols_best.d = ordered_resols[par_depth].pop();
			}
		}
		//printf("%i equivalent are waiting, putting em back!\n", ordered_resols_buffer.getSize());
		for(i=0;i< ordered_resols_buffer.getSize();i++)	ordered_resols[par_depth].insert(ordered_resols_buffer[i]);
		ordered_resols_buffer.toMemfree();



		curalias = U32_Alias(regular.data.deref_alias(regular.data.find(ordered_resols_best.d.d[0])));

		//printf("would do %X\n", ordered_resols_best.d.k);
		if ((i = ordered_resols_best.d.k & 7) != 0){
			regular.insertParentNodeAt(symbufcur++, (uint32_t)'i', curalias , i);
			i = symbufcur-1;
		}else i = ordered_resols_best.d.d[0];
		if (auto ite = regular.mkIterator(i)){
			switch(ordered_resols_best.d.k){
				case 0x133:{ // Addition
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, ((*fite) & 7) | LFHSCRIPT_ADD_8BITS, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x143:{ // Substraction
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, ((*fite) & 7) | LFHSCRIPT_SUB_8BITS, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x113:{ // Multiplication
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, ((*fite) & 7) | LFHSCRIPT_MULT_8BITS, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x123:{ // Division
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, ((*fite) & 7) | LFHSCRIPT_DIVI_8BITS, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x303:{ // <
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_LESSTHAN, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x313:{ // >
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_GRTRTHAN, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x323:{ // <=
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_LESSEQTHAN, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x333:{ // >=
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_GRTREQTHAN, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x343:{ // ==
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_EQUAL, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x353:{ // !=
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = 0x0104;
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_CHECK_UNEQUAL, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x403:{ // &
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_AND, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x423:{ // |
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (*fite);
				if (mkCommand_matchtypes_routine(fite, LFHSCRIPT_OR, strace.getSize())) {printf("constant types mismatch\n");exit(0);}
			}break; case 0x433:{ // && AND
				auto fite = ite.mkPreordernext();
				if (fite() > strace.getSize()) *fite =0;
				fite++;	*fite = LFHSCRIPT_GOTO_IFFALSE;
				fite++;	if (fite() > strace.getSize()) *fite =0;
				*ite = 0x0104;
			}break; case 0x443:{ // || OR
				auto fite = ite.mkPreordernext();
				if (fite() > strace.getSize()) *fite =0;
				fite++;	*fite = LFHSCRIPT_GOTO_IFTRUE;
				fite++;	if (fite() > strace.getSize()) *fite =0;
				*ite = 0x0104;
			}break; case 0x505:{ // OMG ternary expression!
				auto fite = ite.mkPreordernext();
				if (fite() > strace.getSize()) *fite =0;
				fite++; *fite = LFHSCRIPT_GOTO_IFTRUE_TERNARY;
				fite++; *ite = *fite;
				if (fite() > strace.getSize()) *fite =0;
				fite++; *fite = LFHSCRIPT_GOTO_IFFALSE_TERNARY;
				fite++; if (*ite != *fite) {errorSet("Ternary expression has different types for true and false expressions!\n"); return(fout);}
				if (fite() > strace.getSize()) *fite =0;
			}break; case 0x903:{ // e,e
				auto fite = ite.mkPreordernext(); fite.moveNextForward();
				(*ite) = (uint16_t)'c';
			}break; case 0x900:{ // c,e
				//regular.show();
				ite.removeNext();
				//regular.show();
				//ite.show();
				ite.moveNextUnder(-2);
				//regular.show();
				//ite.show();
			}break; case 0x103:{ // silly paranteris (e)
				auto fite = ite.mkPreordernext();
				*fite = 0;
				fite++;
				(*ite) = (*fite);
				if (fite() > strace.getSize()) *fite = 0;
				fite++;*fite = 0;
				par_depth--; // mark parentesis undoing
			}break; case 0x3:{ // function f(void)
				auto fite = ite.mkPreordernext();

				par_depth--; // mark parentesis undoing


				if (auto fncite = fnccodes.mkIteratorFromFirstAlias(constants[ fite() ].int32_val)){
					while(fncite->getSize() != 1){
						if (!(fncite++)) {
							errorNoArgumentMatch_routine(fncite(),"with no argument");
							return fout;
						}
					}
					*fite = 0;fite++;
					*fite = 0;fite++;
					*fite = fncite.getAlias();
					*ite = (*fncite)[0];
				}else{
					fnccodes.show();
					printf("thithithit is imppposible! alias %X\n", constants[ fite() ].int32_val);
					uint32_t mmm = fnccodes.find(U32_Alias(constants[ fite() ].int32_val));
					printf("thithithit is imppposible! ptroff %i\n", mmm);
					printf("thithithit is imppposible! %s\n", fnccodes.deref_key(mmm).c_str());
					exit(1);
				}
			}break; case 0x4:{ // 1 arg function  f(e)
				par_depth--; // mark parentesis undoing
				auto fite = ite.mkPreordernext();
				*fite = 0;
				uint32_t fpos = fite();
				fite++; *fite = 0;
				fite++;
				if (auto fncite = fnccodes.mkIteratorFromFirstAlias(constants[ fpos ].int32_val)){
					while((fncite->getSize() != 2)||((*fncite)[1] != *fite)){
						if (!(fncite++)) {

							errorNoArgumentMatch_routine(fncite(),"with one argument");
							return fout;
						}
					}
					if (fite() > strace.getSize()) *fite = 0;
					fite++;
					*fite = fncite.getAlias();
					*ite = (*fncite)[0];
				}else{
					fnccodes.show();
					printf("thithithit is imppposible! alias %i:%X\n", fpos, constants[ fpos ].int32_val);
					uint32_t mmm = fnccodes.find(U32_Alias(constants[ fpos ].int32_val));
					printf("thithithit is imppposible! ptroff %i\n", mmm);
					printf("thithithit is imppposible! %s\n", fnccodes.deref_key(mmm).c_str());
					exit(1);
				}
			}break;case 0x14:{ // 2+ arg function  f(e,e,e,e...)
				par_depth--; // mark parentesis undoing
				auto fite = ite.mkPreordernext();
				*fite = 0;
				uint32_t fpos = fite();
				fite++; *fite = 0;
				fite++; *fite = 0;
				j=0;
				if (auto site = fite.mkPreordernext()) do {
					if (*site == ',') break;
					j++;
				}while(site++);
				if (auto fncite = fnccodes.mkIteratorFromFirstAlias(constants[ fpos ].int32_val)){
					while(true){
						if (fncite->getSize() == j+1){
								j=1;
							if (auto site = fite.mkPreordernext()) do {
								if (*site == ',') {j=0; break;}
								if ((*fncite)[j++] != *site) break;
							}while(site++);
							if (j == 0) break;
						}
						if (!(fncite++)) {
							errorNoArgumentMatch_routine(fncite(),"with 2+ argument");
							return fout;
						}
					}

					if (auto site = fite.mkPreordernext()) do {
						if (*site == ',') {*site = 0; break;}
						if (site() > strace.getSize()) *site = 0;
					}while(site++);

					fite++; *fite = 0; fite++;
					*fite = fncite.getAlias();
					*ite = (*fncite)[0];

				}else{
					fnccodes.show();
					printf("thithithit is imppposible! alias %i:%X\n", fpos, constants[ fpos ].int32_val);
					uint32_t mmm = fnccodes.find(U32_Alias(constants[ fpos ].int32_val));
					printf("thithithit is imppposible! ptroff %i\n", mmm);
					printf("thithithit is imppposible! %s\n", fnccodes.deref_key(mmm).c_str());
					exit(1);
				}
			}break;default: printf("unregnosized! %X\n", ordered_resols_best.d.k);  exit(1);
			}
		}

		//regular.show();
		// find 7 before insert;
		if (auto ite = regular.mkIterator(0)) {
			if (auto fite = regular.mkIterator(0)) {
				// try to fast forward to 7 before current insertion
				for(partsum =0;partsum<7;partsum++) {
					if (fite() == symbufcur-1) break;
					fite.next();
				}
				if (partsum == 7) {
					while(fite() != symbufcur-1){fite.next();ite.next();}
				}
			}
			i=0;
			partsumbuf[i] =0;
			do{
				partsumbuf[i+1] = partsumbuf[i] + ite();
				for(unsigned int pcur =0; pcur < partial.getSize();pcur++){
					ordered_resols_elem.k = 0;
					//printf("partial code %X\n", partial[pcur] | (((*ite) <= 0x100) ? (*ite) : 0x100) );
					switch( partial[pcur] | (((*ite) <= 0x100) ? (*ite) : 0x100)){
						case 0x00100000 | '!': partial[pcur] = 0x02210000;   // unnary !
					break; case 0x00100000 | '-': partial[pcur] = 0x022D0000;   // unnary -
					break; case 0x00100000 | '*': partial[pcur] = 0x022A0000;   // unnary *
					break; case 0x00100000 | '&': partial[pcur] = 0x02260000;   // unnary &
					break; case 0x00030100: partial.pop_swap(pcur--); // -i detected
					break; case 0x00050100: partial.pop_swap(pcur--); // !i detected
					break; case 0x00040000 | '(': partial[pcur] = 0x01040000; // f( detected
					break; case 0x01040100: partial[pcur] = 0x02040000; // f(i detected
					break; case 0x01040000 | 'c': partial[pcur] = 0x03040000; // f(c detected
					break; case 0x00020100: partial[pcur] = 0x01020000; // (i
					break; case 0x00060000 | ',': partial[pcur] = 0x01C10100; // c,

					break; case 0x00010000 | '?': partial[pcur] = 0x003F0000; // e?
					break; case 0x00010000 | '+': partial[pcur] = 0x002B0000;
					break; case 0x00010000 | '-': partial[pcur] = 0x002D0000;
					break; case 0x00010000 | '*': partial[pcur] = 0x002A0000;
					break; case 0x00010000 | '/': partial[pcur] = 0x002F0000;
					break; case 0x00010000 | '%': partial[pcur] = 0x00250000;

					break; case 0x00010000 | '<': partial[pcur] = 0x003C0000;
					break; case 0x00010000 | '=': partial[pcur] = 0x003D0000;
					break; case 0x00010000 | '>': partial[pcur] = 0x003E0000;
					break; case 0x00010080 | '<': partial[pcur] = 0x013C0000;
					break; case 0x00010080 | '=': partial[pcur] = 0x013D0000;
					break; case 0x00010080 | '>': partial[pcur] = 0x013E0000;
					break; case 0x00010000 | '&': partial[pcur] = 0x00260000;
					break; case 0x00010080 | '&': partial[pcur] = 0x01260000;
					break; case 0x00010000 | '|': partial[pcur] = 0x007C0000;
					break; case 0x00010080 | '|': partial[pcur] = 0x017C0000;
					break; case 0x00010080 | '!': partial[pcur] = 0x01210000;
					break; case 0x00010000 | ',': partial[pcur] = 0x00C10100;


					//break; case 0x00020000 | '+': partial[pcur] = 0x01020000;
					//break; case 0x00020000 | '-': partial[pcur] = 0x01020000;
					//break; case 0x00020000 | '*': partial[pcur] = 0x02020000;
					//break; case 0x00020000 | '/': partial[pcur] = 0x02020000;

					break; case 0x00050000 | 'e': partial[pcur] = 0x01050000;

					break; case 0x003F0100: partial[pcur] = 0x013F0000; // e?e
					break; case 0x013F0000 | ':': partial[pcur] = 0x023F0000; // e?e:
					break; case 0x023F0100: partial[pcur] = ordered_resols_elem.k = 0x505; // e?e:e

					break; case 0x01040000 | ')': ordered_resols_elem.k = 0x3; // f() detected
					break; case 0x02040000 | ')': ordered_resols_elem.k = 0x4; // f(i) detected
					break; case 0x03040000 | ')': ordered_resols_elem.k = 0x14; // f(c) detected
					break; case 0x01020000 | ')': ordered_resols_elem.k = 0x103; // (e) detected

					break; case 0x002A0100: ordered_resols_elem.k = 0x113;// e * e detected
					break; case 0x002F0100: ordered_resols_elem.k = 0x123;// e / e detected
					break; case 0x002B0100: ordered_resols_elem.k = 0x133;// e + e detected
					break; case 0x002D0100: ordered_resols_elem.k = 0x143;// e - e detected
					break; case 0x00250100: ordered_resols_elem.k = 0x153;// e % e detected

					break; case 0x02210100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary ! e detected
					break; case 0x022D0100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;  // unnary - e detected
					break; case 0x022A0100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary * e detected
					break; case 0x02260100: partial.pop_swap(pcur--);//ordered_resols_elem.k = 0x303;   // unnary & e detected

					break; case 0x003C0100: ordered_resols_elem.k = 0x303; // e < e detected
					break; case 0x013C0100: ordered_resols_elem.k = 0x323; // e <= e detected
					break; case 0x003D0100: ordered_resols_elem.k = 0x733; // e = e detected
					break; case 0x013D0100: ordered_resols_elem.k = 0x343; // e == e detected
					break; case 0x01210000: ordered_resols_elem.k = 0x353; // e != e detected
					break; case 0x003E0100: ordered_resols_elem.k = 0x313; // e > e detected
					break; case 0x013E0100: ordered_resols_elem.k = 0x333; // e >= e detected
					break; case 0x00260100: ordered_resols_elem.k = 0x403; // e & e detected
					break; case 0x01260100: ordered_resols_elem.k = 0x433; // e && e detected
					break; case 0x007C0100: ordered_resols_elem.k = 0x423; // e | e detected
					break; case 0x017C0100: ordered_resols_elem.k = 0x443; // e || e detected

					break; case 0x00C10100: ordered_resols_elem.k = 0x903; // e,e detected
					break; case 0x01C10100: ordered_resols_elem.k = 0x900; // c,e detected
					//break; case 0x02020000 | 'f': ordered_resols_elem.k = 0x53; // e * e detected
					break; default: partial.pop_swap(pcur--);
					}
					if (ordered_resols_elem.k != 0){
						//printf("complete... but? %i %i %i\n", i- ((ordered_resols_elem.k & 7) -1), partsum, i);
						if ((j = (ordered_resols_elem.k & 7)) == 0){ // ah... special case
							j = (ordered_resols_elem.k == 0x900) ? 3 : 4;
						}
						if ((i >= partsum)&&(i- (j -1) <= partsum ))  { // now must contain new position!
							partial.pop_swap(pcur--);
							ordered_resols_elem.d[0] = partsumbuf[i - j+2] - partsumbuf[i- j+1];
							ordered_resols_elem.d[1] = ite();
							ordered_resols_elem.d[2] = partsumbuf[i+1] - partsumbuf[i- j+1];
							//printf("newcomplete!:\n"); ordered_resols_elem.show();
							ordered_resols[par_depth].insert(ordered_resols_elem);
						}
					}
				}
				//printf("S[%i] = %c: %i partial still satisfied\n", symbufcur, (char) (*ite), partial.getSize());
				// add rules resolutions

				if ((*ite) >= 0x100) partial.push_back(0x10000); // implicit conversion to expression
				else switch(*ite){
					case '(': partial.push_back(0x20000); partial.push_back(0x100000);
				break; case 'F': partial.push_back(0x40000);
				break; case 'c': partial.push_back(0x60000);
				break; default:  partial.push_back(0x100000);
				}
				//printf("at pos %i: %i partial, %i complete\n", symbufcur, partial.getSize(), ordered_resols[par_depth].getSize());
				if (i++ == 14) break;
			}while(ite.next());
			partial.toMemfree();
		}
		par_depth = ordered_resols.getSize()-1;
	}

	// check if return value is found

	//regular.show();

	Vector<uint16_t> scriptbuf;
	// insert constant insertion commands
	for(uint32_t i = 0; i < strace.getSize();i++){
		if (auto ite = regular.mkIterator(i+1)){
			if (((*ite) >= 0x100)&&((*ite) < 0xF000)) { // that's not great...
				(*ite) = (((*ite) & 0x002) != 0) ? ((((*ite) & 0x001) != 0) ? LFHSCRIPT_LOAD_CONST_64BITS : LFHSCRIPT_LOAD_CONST_32BITS ) : LFHSCRIPT_LOAD_CONST_16BITS;
			}
		}
	}

	if (auto ite = regular.mkIterator(0)) {
		//printf("%i:%c%X\t", ite(), (*ite) &  255, (*ite) >> 8 );
		if (ite.next()){
			if (ite() >= strace.getSize()) *ite = 0;
			if (ite.next()){
				errorSet("parser failed!\n"); return(fout);
			}
		}else{
			errorSet("parser failed!\n"); return(fout);
		}
	}

	//constants.show();
	symbufcur = 0;
	myHashmap<uint32_t, KeyElem< Vector<uint32_t> , Vector<uint32_t> > > addrtofill_purgatory;
	Vector<uint32_t> addr_orfill; // aka goto if true
	Vector<uint32_t> addr_andfill; // aka goto if false
	bool wasboolean;
	int startaddr;
	if (auto ite = regular.mkIterator(0)) while(ite.postordernext()){
		uint32_t i;
		wasboolean= false;
		startaddr = scriptbuf.getSize();
		switch(*ite){
		       case 0: wasboolean = true; // nothing done, look for next actual operation
		break; case LFHSCRIPT_LOAD_CONST_16BITS: scriptbuf.push_back(*ite); scriptbuf.push_back(constants[ ite() ].int16_val); symbufcur++;
		break; case LFHSCRIPT_LOAD_CONST_32BITS: scriptbuf.upSize(scriptbuf.getSize()+3); scriptbuf.last(3) = (*ite); *(uint32_t*)(&scriptbuf.last(2)) = constants[ite()].int32_val; symbufcur++;
		break; case LFHSCRIPT_LOAD_CONST_64BITS: scriptbuf.upSize(scriptbuf.getSize()+5); scriptbuf.last(5) = (*ite); *(uint64_t*)(&scriptbuf.last(4)) = constants[ite()].int64_val; symbufcur++;
		break; case LFHSCRIPT_GOTO_IFTRUE: wasboolean = true;
			addr_orfill.push_back(scriptbuf.getSize()); scriptbuf.push_back(*ite);
			addrtofill_purgatory[ite.getParentAlias()].k.toMemmove(addr_orfill); scriptbuf.push_back(0);
			for(i=0; i < addr_andfill.getSize();i++){
				scriptbuf[addr_andfill[i]] = LFHSCRIPT_GOTO_IFFALSE_TERNARY;
				scriptbuf[addr_andfill[i]+1] = scriptbuf.getSize();
			}
			addr_andfill.toMemfree();
		break; case LFHSCRIPT_GOTO_IFFALSE: wasboolean = true;
			addr_andfill.push_back(scriptbuf.getSize()); scriptbuf.push_back(*ite);
			addrtofill_purgatory[ite.getParentAlias()].d.toMemmove(addr_andfill); scriptbuf.push_back(0);
			for(i=0; i < addr_orfill.getSize();i++){
				scriptbuf[addr_orfill[i]] = LFHSCRIPT_GOTO_IFTRUE_TERNARY;
				scriptbuf[addr_orfill[i]+1] = scriptbuf.getSize();
			}
			addr_orfill.toMemfree();
		break; case LFHSCRIPT_GOTO_IFTRUE_TERNARY: // ?  TODO scriptbuf.push_back(*ite); addrtofill.push_back(scriptbuf.getSize()); scriptbuf.push_back(0);
			addr_andfill.push_back(scriptbuf.getSize()); scriptbuf.push_back(*ite);
			addrtofill_purgatory[ite.getParentAlias()].d.toMemmove(addr_andfill); scriptbuf.push_back(0);
			for(i=0; i < addr_orfill.getSize();i++){
				scriptbuf[addr_orfill[i]] = LFHSCRIPT_GOTO_IFTRUE_TERNARY;
				scriptbuf[addr_orfill[i]+1] = scriptbuf.getSize();
			}
			addr_orfill.toMemfree();
		break; case LFHSCRIPT_GOTO_IFFALSE_TERNARY: // : TODO scriptbuf.push_back(*ite); addrtofill.push_back(scriptbuf.getSize()); scriptbuf.push_back(0);
			addr_andfill.push_back(scriptbuf.getSize()); scriptbuf.push_back(LFHSCRIPT_GOTO);
			addrtofill_purgatory[ite.getParentAlias()].k.toMemmove(addr_andfill); scriptbuf.push_back(0);
			for(i=0; i < addrtofill_purgatory[ite.getParentAlias()].d.getSize();i++){
				scriptbuf[addrtofill_purgatory[ite.getParentAlias()].d[i]] = LFHSCRIPT_GOTO_IFFALSE_TERNARY;
				scriptbuf[addrtofill_purgatory[ite.getParentAlias()].d[i]+1] = scriptbuf.getSize();
			}
			addrtofill_purgatory[ite.getParentAlias()].d.toMemfree();
		break; default:
			if ((*ite) < 0xF000) if (fnccodes[ U32_Alias(*ite)].getSize() == 1) symbufcur++;
			scriptbuf.push_back(*ite);
		break;}
		if (!wasboolean){
			for(i=0; i < addr_orfill.getSize();i++) scriptbuf[addr_orfill[i]+1] = startaddr;
			for(i=0; i < addr_andfill.getSize();i++) scriptbuf[addr_andfill[i]+1] = startaddr;
			addr_orfill.toMemfree();
			addr_andfill.toMemfree();
		}//((3>2)||(boolgen(5)))&&(boolgen(3))

		//if (addrtofill_purgatory.getSize() != 0) {printf("%i (aka 0x%X) to check in:\n", ite.getAlias().value, ite.getAlias().value);
		//	addrtofill_purgatory.show();
		//}
		if ((i = addrtofill_purgatory.find(ite.getAlias().value)) != 0xFFFFFFFF){
			addr_orfill.toMemappend(addrtofill_purgatory.deref(i).k);
			addr_andfill.toMemappend(addrtofill_purgatory.deref(i).d);
			addrtofill_purgatory.erase_from_iterator(i);
		}

	}
	for(i=0; i < addr_orfill.getSize();i++) scriptbuf[addr_orfill[i]+1] = scriptbuf.getSize();
	for(i=0; i < addr_andfill.getSize();i++) scriptbuf[addr_andfill[i]+1] = scriptbuf.getSize();
	addr_orfill.toMemfree();
	addr_andfill.toMemfree();

	fout.code.setSize(scriptbuf.getSize()+1);
	fout.code[0] = symbufcur;
	for(unsigned int i =0 ; i< scriptbuf.getSize(); i++) fout.code[i+1] = scriptbuf[i];
return fout;}

void ScriptCompiler::loadDeclarations(const Vector< KeyElem< KeyElem<string, uint32_t>, Tuple<uint32_t> > >& dacst){
	for(unsigned int i =0 ; i< dacst.getSize();i++) fnccodes.createEntryWithAlias(dacst[i].k.k, dacst[i].k.d) = dacst[i].d;
}

void ScriptCompiler::saveDeclarations(Vector< KeyElem< KeyElem<string, uint32_t> , Tuple<uint32_t> > >& dacst) const{
	if (auto ite = fnccodes.mkIterator()) do{
		dacst.push_back();
		dacst.last().k.k = ite();
		dacst.last().k.d = ite.getAlias();
		dacst.last().d = (*ite);
	}while(ite++);
}

Ambiguous ScriptCompiler::executeCode(const uint16_t* code, uint32_t length, ScriptScope* fncscopr){
	Tuple<Ambiguous> scope; scope.setSize(code[0]);
	Ambiguous* ptr = scope();
	char* error;
	for(unsigned int i =1; i < length; i++){
		//for(Ambiguous* sptr = scope.data; sptr != ptr; sptr++) printf("%X\t", sptr->uint32_val);
		//printf("(scope at %i)\n", i);
		//printf("torun code %X\n", code[i]);
		if (code[i] <= (uint32_t) LFHSCRIPT_LAST_OPERATION) { error = externfnc[code[i]](ptr);
			if (error != NULL){
				printf("Error %s\n",error );
				delete[](error);
				return Ambiguous();
			}
		}else switch(code[i] ){
		//	case LFHSCRIPT_CHECK_EQUAL_8BITS: ptr--; ptr[-1].int32_val = ( (ptr[-1].uint8_val ^ ptr->uint8_val) != 0) ? 1 : 0;
		//break; case LFHSCRIPT_CHECK_EQUAL_16BITS: ptr--; ptr[-1].int32_val = ( (ptr[-1].uint16_val ^ ptr->uint16_val) != 0) ? 1 : 0;
		//break; case LFHSCRIPT_CHECK_EQUAL_32BITS: ptr--; ptr[-1].int32_val = ( (ptr[-1].uint32_val ^ ptr->uint32_val) != 0) ? 1 : 0;
		//break; case LFHSCRIPT_CHECK_EQUAL_64BITS: ptr--; ptr[-1].int32_val = ( (ptr[-1].uint64_val ^ ptr->uint64_val) != 0) ? 1 : 0;

			case LFHSCRIPT_CHECK_EQUAL: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val == ptr->int64_val) ? 1 : 0;
		break; case LFHSCRIPT_CHECK_UNEQUAL: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val != ptr->int64_val) ? 1 : 0;
		break; case LFHSCRIPT_CHECK_LESSTHAN: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val < ptr->int64_val) ? 1 : 0;
		break; case LFHSCRIPT_CHECK_GRTRTHAN: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val > ptr->int64_val) ? 1 : 0;
		break; case LFHSCRIPT_CHECK_LESSEQTHAN: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val <= ptr->int64_val) ? 1 : 0;
		break; case LFHSCRIPT_CHECK_GRTREQTHAN: ptr--; ptr[-1].int64_val = (ptr[-1].int64_val <= ptr->int64_val) ? 1 : 0;

		break; case LFHSCRIPT_64BITS_TO_F32BITS: ptr[-1].float_val = (float) ptr[-1].uint64_val;
		break; case LFHSCRIPT_64BITS_TO_F64BITS: ptr[-1].double_val = (double) ptr[-1].uint64_val;
		break; case LFHSCRIPT_F64BITS_TO_F32BITS: ptr[-1].float_val = (float) ptr[-1].double_val;
		break; case LFHSCRIPT_F32BITS_TO_F64BITS: ptr[-1].double_val = (double) ptr[-1].float_val;

		break; case LFHSCRIPT_LOAD_CONST_ZERO: ptr->uint64_val = 0; ptr++;
		break; case LFHSCRIPT_LOAD_CONST_16BITS: ptr->uint64_val = 0; ptr->uint16_val = code[++i]; ptr++;
		break; case LFHSCRIPT_LOAD_CONST_32BITS: ptr->uint64_val = 0; ptr->uint32_val = code[++i]; ptr->uint32_val |= (((uint32_t)code[++i])<<16); ptr++;
		break; case LFHSCRIPT_LOAD_CONST_64BITS: ptr->uint64_val = 0; ptr->uint64_val = code[++i]; ptr->uint32_val |= (((uint32_t)code[++i])<<16); ptr->uint64_val |= (((uint64_t)code[++i])<<32); ptr->uint64_val |= (((uint64_t)code[++i])<<48); ptr++;

		break; case LFHSCRIPT_ADD_8BITS: ptr--; ptr[-1].int8_val += ptr->uint8_val;
		break; case LFHSCRIPT_ADD_16BITS: ptr--; ptr[-1].int16_val += ptr->uint16_val;
		break; case LFHSCRIPT_ADD_32BITS: ptr--; ptr[-1].int32_val += ptr->uint32_val;
		break; case LFHSCRIPT_ADD_64BITS: ptr--; ptr[-1].int64_val += ptr->uint64_val;
		break; case LFHSCRIPT_OR: ptr--; ptr[-1].int64_val |= ptr->int64_val;
		break; case LFHSCRIPT_ADD_F32BITS: ptr--; ptr[-1].float_val += ptr->float_val;
		break; case LFHSCRIPT_ADD_F64BITS: ptr--; ptr[-1].double_val += ptr->double_val;

		break; case LFHSCRIPT_SUB_8BITS: ptr--; ptr[-1].int8_val -= ptr->uint8_val;
		break; case LFHSCRIPT_SUB_16BITS: ptr--; ptr[-1].int16_val -= ptr->uint16_val;
		break; case LFHSCRIPT_SUB_32BITS: ptr--; ptr[-1].int32_val -= ptr->uint32_val;
		break; case LFHSCRIPT_SUB_64BITS: ptr--; ptr[-1].int64_val -= ptr->uint64_val;
		break; case LFHSCRIPT_XOR: ptr--; ptr[-1].int64_val ^= ptr->int64_val;
		break; case LFHSCRIPT_SUB_F32BITS: ptr--; ptr[-1].float_val -= ptr->float_val;
		break; case LFHSCRIPT_SUB_F64BITS: ptr--; ptr[-1].double_val -= ptr->double_val;


		break; case LFHSCRIPT_MULT_8BITS: ptr--; ptr[-1].int8_val *= ptr->uint8_val;
		break; case LFHSCRIPT_MULT_16BITS: ptr--; ptr[-1].int16_val *= ptr->uint16_val;
		break; case LFHSCRIPT_MULT_32BITS: ptr--; ptr[-1].int32_val *= ptr->uint32_val;
		break; case LFHSCRIPT_MULT_64BITS: ptr--; ptr[-1].int64_val *= ptr->uint64_val;
		break; case LFHSCRIPT_AND: ptr--; ptr[-1].int64_val &= ptr->int64_val;
		break; case LFHSCRIPT_MULT_F32BITS: ptr--; ptr[-1].float_val *= ptr->float_val;
		break; case LFHSCRIPT_MULT_F64BITS: ptr--; ptr[-1].double_val *= ptr->double_val;

		break; case LFHSCRIPT_DIVI_8BITS: ptr--; ptr[-1].int8_val /= ptr->uint8_val;
		break; case LFHSCRIPT_DIVI_16BITS: ptr--; ptr[-1].int16_val /= ptr->uint16_val;
		break; case LFHSCRIPT_DIVI_32BITS: ptr--; ptr[-1].int32_val /= ptr->uint32_val;
		break; case LFHSCRIPT_DIVI_64BITS: ptr--; ptr[-1].int64_val /= ptr->uint64_val;
		break; case LFHSCRIPT_DIVI_F32BITS: ptr--; ptr[-1].float_val /= ptr->float_val;
		break; case LFHSCRIPT_DIVI_F64BITS: ptr--; ptr[-1].double_val /= ptr->double_val;

		break; case LFHSCRIPT_GOTO: i++; i = code[i];
		break; case LFHSCRIPT_GOTO_IFFALSE_TERNARY_FINAL: i++; if (ptr[-1].uint32_val == 0) i = code[i]; ptr--;
		break; case LFHSCRIPT_GOTO_IFTRUE: i++; if (ptr[-1].uint32_val != 0) i = code[i]; else ptr--;
		break; case LFHSCRIPT_GOTO_IFFALSE: i++; if (ptr[-1].uint32_val == 0) i = code[i]; else  ptr--;
		break; case LFHSCRIPT_GOTO_IFTRUE_TERNARY: i++; if (ptr[-1].uint32_val != 0) i = code[i]; ptr--;
		break; case LFHSCRIPT_GOTO_IFFALSE_TERNARY: i++; if (ptr[-1].uint32_val == 0) i = code[i];  ptr--;
		}
	}
return scope[0];}



const ScriptCompiler& ScriptCompiler::show(FILE* f, int level) const{
	fprintf(f,"Script scope, list of functions:");
	if (auto ite = fnccodes.mkIterator()) do{
		fprintf(f,"\n%i", ite.getAlias());
		fprintf(f,":%X", (*ite)[0]);
		fprintf(f,"\t%s",  ite().c_str());
		for(unsigned int i =1; i < (*ite).getSize();i++) fprintf(f,"\t%X", (*ite)[i]);
	}while(ite++);
	fprintf(f,"\n");
return *this;}

Ambiguous Script::execute() const{  // std::function< void(uint32_t, Ambiguous*&) > externfnc
	if (code.getSize() == 0){ // assuming this just got compiled with unmanaged error handling, raise error and exit
		//if ((funchost.errbuffer) != NULL) printf("%s", funchost.errbuffer);
		//delete[](funchost.errbuffer); funchost.errbuffer = NULL;
		return Ambiguous();
	}
return funchost.executeCode(code(), code.getSize());}



const Script& Script::show(FILE* f, int level) const{
	fprintf(f,"Script with length %i:\n", code.getSize());
	for(uint32_t i = 0; i < code.getSize(); i++) {
		switch(code[i] & 0xFFF8){
			case LFHSCRIPT_ADD_8BITS:  fprintf(f,"%i:ADD%X\t", i,code[i] & 7);
		break; case LFHSCRIPT_SUB_8BITS:  fprintf(f,"%i:SUB%X\t", i, code[i] & 7);
		break; case LFHSCRIPT_MULT_8BITS:  fprintf(f,"%i:MULT%X\t", i, code[i] & 7);
		break; case LFHSCRIPT_DIVI_8BITS:  fprintf(f,"%i:DIVI%X\t", i, code[i] & 7);
		break; case LFHSCRIPT_LOAD_CONST_ZERO:
			if (code[i] & 4){
				fprintf(f,(code[i] & 1) ? "%i:(double)\t" : "%i:(float)\t",i);
			}else{
				if (code[i] & 2){
					if (code[i++] & 1){
						fprintf(f,"%i:LOAD %i\t",i-1, (uint32_t) *(uint64_t*)&code[i]);
						i+=3;
					}else {fprintf(f,"%i:LOAD %i\t", i-1,*(uint32_t*)&code[i]); i++;}
				}else if (code[i] & 1){
					fprintf(f,"%i:LOAD %i\t", i, code[i+1]);
				}else{
					fprintf(f,"%i:LOADZERO\t",i);
				}
			}
		break; case LFHSCRIPT_GOTO_:
			switch(code[i]){
				case LFHSCRIPT_GOTO: fprintf(f,"%i:GOTO %i\t",i, code[i+1] +1);
			break; case LFHSCRIPT_GOTO_IFFALSE_TERNARY_FINAL: fprintf(f,"%i:GOTO_IF_F_AC %i\t",i, code[i+1] +1);
			break; case LFHSCRIPT_GOTO_IFTRUE: fprintf(f,"%i:GOTO_IF_T %i\t",i, code[i+1] +1);
			break; case LFHSCRIPT_GOTO_IFFALSE: fprintf(f,"%i:GOTO_IF_F %i\t",i, code[i+1] +1);
			break; case LFHSCRIPT_GOTO_IFTRUE_TERNARY: fprintf(f,"%i:GOTO_IF_T_C %i\t",i, code[i+1] +1);
			break; case LFHSCRIPT_GOTO_IFFALSE_TERNARY: fprintf(f,"%i:GOTO_IF_F_C %i\t",i, code[i+1] +1);
			}
			i++;
		break; default:
			if (code[i] <= 0xF000) {
				uint32_t daite = funchost.fnccodes.find(U32_Alias(code[i]));
				if (daite == 0xFFFFFFFF) fprintf(f,"%i:UNKNOWN%X\t",i, code[i]);
				else fprintf(f,"%i:%s\t",i, funchost.fnccodes.deref_key(daite).c_str());
			}else fprintf(f,"%u:FUNC%X\t",i, code[i]);
		}
	}
	fprintf(f,"\n");
return *this;}

void ScriptScopeOld::initdictionary(){
	dico[string("int")] = 0;
	dico[string("float")] = 0;
	dico[string("return")] = 0;
	dico[string("if")] = 0;
	dico[string("else")] = 0;
	}
void ScriptScopeOld::finitdictionary(){ExOp::toMemfree(dico);}
void ScriptScopeOld::compileScript(FILE* f){
    int stacksize=4;
    myHashmap<string, unsigned int> variables;
    Vector<char> tmp_code;
    char buffer[65536];
    uint32_t i;
    uint32_t lineNo =1;
    char* cur;
    Vector<uint32_t> parse_state;
    parse_state.push_back(0);
    while(feof(f) == false){
        if (1 == fscanf(f, "%[a-zA-Z0-9_] \t\n",buffer)){
            i = dico.find(string(buffer));
            if (i == 0xFFFFFFFF) {
                i = variables.find(buffer);
                if (i == 0xFFFFFFFF) {
                    switch(parse_state.last()){
                    case 1:
                        if ((buffer[0] >= '0')&&(buffer[0] <= '9')) {LFH_exit("Line %i: Variable name cannot start with a number\n", lineNo);}
                        printf("successfully declared int %s\n", buffer);
                        variables[buffer] = stacksize;
                        stacksize+= 4;
                    break;
                    case 2:
                        if ((buffer[0] >= '0')&&(buffer[0] <= '9')) {LFH_exit("Line %i: Variable name cannot start with a number\n", lineNo);}
                        printf("successfully declared float %s\n", buffer);
                        variables[buffer] = stacksize | 0x80000000;
                        stacksize+= 4;
                    break;
                    case 3:

                    break;
                    default:
                    printf("unknown token: %s\n", buffer);
                    }
                }else{
                    switch(parse_state.last()){
                        case 1:
                        case 2:
                            LFH_exit("Line %i: variable '%s' is already declared!\n", lineNo, buffer);
                        break;
                        default:
                            printf("Using variable %s\n", buffer);
                    }
                }
            }else if (dico.deref(i) == 0) {
                switch(buffer[0]){
                    case 'e':
                        if (strcmp(buffer+1,"lse")==0){
                            printf("else token!\n");
                        }
                    break;
                    case 'f':
                        if (strcmp(buffer+1,"loat")==0){
                            if (parse_state.last() != 0) {LFH_exit("Line %i: Unknown token before 'int' token\n", lineNo);}
                            parse_state.last() = 2;
                        }
                    break;
                    case 'i':
                        if (strcmp(buffer+1,"f")==0) {
                            printf("if token!\n");
                        }else if (strcmp(buffer+1,"nt")==0) {
                            if (parse_state.last() != 0) {LFH_exit("Line %i: Unknown token before 'int' token\n", lineNo);}
                            parse_state.last() = 1;
                        }
                    break;
                    case 'r':
                        if (strcmp(buffer+1,"eturn")==0){
                            printf("return token!\n");
                        }
                    break;
                    }
            }else{
                printf("dico token: %s\n", buffer);
            }
        }else if (1 == fscanf(f, "%[!-/:-?[-]{-~]",buffer)) {
            cur = buffer;
            while(*cur != '\0'){
                switch(*cur++){
                case ';': parse_state.last() =0;
                break;

                case '=':
                    if (*cur == '\0') parse_state.last() =3;
                    else {fprintf(stderr,"Line %i: Operator token unrecogninized (%s)\n", lineNo, buffer); LFH_exit(1);}
                break;
                case '[':
                    parse_state.push_back(4);
                break;
                case ']':
                    if (parse_state.last() != 4) {fprintf(stderr,"Line %i: token ']' does not match a '['\n", lineNo); LFH_exit(1);}
                    parse_state.pop_back();
                break;
                case '(':
                    parse_state.push_back(5);
                break;
                case ')':
                    if (parse_state.last() != 5) {fprintf(stderr,"Line %i: token ']' does not match a '['\n", lineNo); LFH_exit(1);}
                    parse_state.pop_back();
                break;
                case '+':
                    switch(*cur){
                        case '=': cur++;printf("operator\n");
                        break;
                        case '+': cur++;printf("operator\n");
                        break;
                        default:
                        printf("operator\n");
                    }
                break;
                case '-':
                    switch(*cur){
                        case '=': cur++;printf("operator\n");
                        break;
                        case '-': cur++;printf("operator\n");
                        break;
                        default:printf("operator\n");
                    }
                break;
                case '<':
                    switch(*cur){
                        case '=': cur++;printf("operator\n");
                        break;
                        case '<': cur++;printf("operator\n");
                        break;
                        default:printf("operator\n");
                    }
                break;
                case '>':
                    switch(*cur){
                        case '=': cur++;printf("operator\n");
                        break;
                        case '>': cur++;printf("operator\n");
                        break;
                        default:printf("operator\n");
                    }
                break;
                default:
                    printf("unknown sep: %s\n", buffer);
                }
            }
        }else {
            i = fgetc(f);
            if (i == '\n') lineNo++;
        }
    }

    stack = new char[stacksize]; LFH_NICE_ALLOCERROR(stack ,"")
    stacksize = tmp_code.getSize();
    script = ( stacksize) ? new char[ stacksize] : NULL; LFH_NICE_ALLOCERROR(script ,"")
    cursor = script;
}
uint32_t ScriptScopeOld::operator()(){
    switch(*cursor++){
        case 0: return 0;
        case 1: cursor+= 4; return *(uint32_t*)(cursor-4); // call external compiled function
        case 2: // copy 4-byte variable
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) = *(uint32_t*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 3: // set 4-byte variable
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) = *(uint32_t*)(cursor+4);
            cursor+= 8;
        break;
        case 4: // add 4-byte int
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) += *(uint32_t*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 5: // sub 4-byte int
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) -= *(uint32_t*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 6: // sub 4-byte int
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) *= *(uint32_t*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 7: // sub 4-byte int
            *(uint32_t*)(stack + (*(uint32_t*)cursor)) /= *(uint32_t*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 8: // add 4-byte int
            *(float*)(stack + (*(uint32_t*)cursor)) += *(float*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 9: // sub 4-byte int
            *(float*)(stack + (*(uint32_t*)cursor)) -= *(float*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 10: // sub 4-byte int
            *(float*)(stack + (*(uint32_t*)cursor)) *= *(float*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 11: // sub 4-byte int
            *(float*)(stack + (*(uint32_t*)cursor)) /= *(float*)(stack + (*(uint32_t*)(cursor+4)));
            cursor+= 8;
        break;
        case 16: // compare 4-byte int
            if (*(int*)(stack + (*(uint32_t*)cursor)) == *(int*)(stack + (*(uint32_t*)(cursor+4)))) compbit =0;
            else compbit = (*(int*)(stack + (*(uint32_t*)cursor)) < *(int*)(stack + (*(uint32_t*)(cursor+4)))) ? 1:2;
        break;
        case 17: // jump
            cursor = script + *(uint32_t*)(stack + (*(uint32_t*)cursor));
        break;
        case 18: // jump if eq
            cursor = ( (compbit &3) == 0) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
        case 19: // jump if neq
            cursor = (compbit &3) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
        case 20: // jump if leq
            cursor = ( (compbit &2) == 0) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
        case 21: // jump if lth
            cursor = ( (compbit &1) == 1) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
        case 22: // jump if leq
            cursor = ( (compbit &1) == 0) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
        case 23: // jump if lth
            cursor = ( (compbit &2) == 2) ? script + *(uint32_t*)(stack + (*(uint32_t*)cursor)) : cursor+4;
        break;
    }
    return 0;
}

/*

    void SQuaternionmk2::operator*=(SQuaternionmk2 other){

    }
    void SQuaternionmk2::show(FILE* f = stdout, int level =0) const{
    switch(level){
        case 0:{
            const float lf = (1.0f / 16384.0f);
        	int rx2 = (int) x;
            int ry2 = (int) y;
            int rz2 = (int) z;
            int nox = (-rx2 * w) >> 15;
            int noy = (-ry2 * w) >> 15;
            int noz = (-rz2 * w) >> 15;
            int nax = (ry2 * rz2) >> 15;
            int nay = (rx2 * rz2) >> 15;
            int naz = (rx2 * ry2) >> 15;
            rx2 = (rx2 * rx2) >> 15;
            ry2 = (ry2 * ry2) >> 15;
            rz2 = (rz2 * rz2) >> 15;
            printf("%.1f\t%.1f\t%.1f\n", lf * (16384 - ry2 -rz2), lf * (-noz +naz), lf * (noy +nay) );
            printf("%.1f\t%.1f\t%.1f\n", lf * (noz +naz), lf * (16384 - rx2 -rz2), lf * (-nox +nax) );
            printf("%.1f\t%.1f\t%.1f\n", lf * (-noy +nay), lf * (-nox +nax), lf * (16384 - rx2 -ry2) );
        }break;
        default:
            fprintf(f,"%i\t%i\t%i\t%i\n",w,x,y,z);
        }
}
*/

/*
	angle::angle(){}
	angle::angle(double value){
		if ((value >= 1.0f)||(value <= 0.0f)) ang =0.0f;
		else ang = (unsigned long long) (value * pow(2.0f, sizeof(unsigned long long) * 8.0f));

		}
	double angle::real(){
			return(cos(((double)ang)*pow(0.5f, sizeof(unsigned long long) * 8.0f) * M_PI*2 ) );
		}
	double angle::imag(){
		  return(sin(((double)ang)*pow(0.5f, sizeof(unsigned long long) * 8.0f) * M_PI*2 ) );
		}*/

/*
static double quadsolve(double* what){
	if (what[2] == 0.0f) {
		if (what[1] == 0.0f) return(what[0] == 0.0f ? 0.0f : -what[0]/ what[1]);
		else return(-what[0]/ what[1]);
	}
	double delta = what[1] * what[1] -4* what[0]* what[2];
	if (delta < 0) return(-100.0f);
	else return( (-what[1] + (what[0] > 0 ? -sqrt(delta) : sqrt(delta))) / (2*what[2]));
}

static double orientedquadsolve(double* what){
	if (what[2] == 0.0f) {
		if (what[1] == 0.0f) return(what[0] == 0.0f ? 0.0f : -what[0]/ what[1]);
		else return(-what[0]/ what[1]);
	}
	double delta = what[1] * what[1] -4* what[0]* what[2];
	if (delta < 0) return(-100.0f);
	else return( (-what[1] - sqrt(delta)) / (2*what[2]));
}

static double negorientedquadsolve(double* what){
	if (what[2] == 0.0f) {
		if (what[1] == 0.0f) return(what[0] == 0.0f ? 0.0f : -what[0]/ what[1]);
		else return(-what[0]/ what[1]);
	}
	double delta = what[1] * what[1] -4* what[0]* what[2];
	if (delta < 0) return(-100.0f);
	else return( (-what[1] + sqrt(delta)) / (2*what[2]));
}*/

	/*
	template < > double PolyThing<double,0>::findMax(double imin, double imax){
		PolyThing<double,0> nd = secondDerivative();
		Vector<complex> zeros = nd.getZeros();

		double best,bestv;
		double tmpv;
		int i;
		best =imin;
		bestv = eval(imin);

		if (imin != imax){
			tmpv = eval(imax);
			if (tmpv > bestv){
				tmpv = bestv;
				best = imax;
			}
		}
		for(i=zeros.size()-1;i>=0;i--){
			if (zeros[i][1] == 0.0f) {
				tmpv = eval(zeros[i][0]);
				if (tmpv > bestv){
					tmpv = bestv;
					best = zeros[i][0];
				}
			}
		}
		return(best);
	}
	template < > Vector<complex,LFHPrimitive::LFHVECTOR_NORMAL> PolyThing<double,0>::getZeros(){
		double tmp, tmp2;
		Vector<complex> _out;
		switch(order){
			case 2:
				tmp = coef[1] / (coef[2] *2);
				tmp2 = tmp*tmp +  coef[0] / coef[2];
				if (tmp2 < 0) {
					tmp2 = sqrt(-tmp2);
					_out.push_back(complex(tmp, tmp2));
					_out.push_back(complex(tmp, -tmp2));
				}else {
					tmp2 = sqrt(tmp2);
					_out.push_back(complex(tmp + tmp2,0.0f));
					_out.push_back(complex(tmp - tmp2,0.0f));
				}
				break;

		}
		return(_out);

	}


*//*
static double systemSolveSph(PolyVector &a, PolyVector &b){
	// norm(a-b) = ar+ br:
	Polything x = a.pos[0] - b.pos[0];
	Polything y = a.pos[1] - b.pos[1];
	Polything z = a.pos[2] - b.pos[2];
	Polything r = a.pos[3] + b.pos[3];

	x.square();
	y.square();
	z.square();
	r.square();

	Polything f = x + y + z - r;
	return(f.solve());
}

static double systemSolveSph(PolyVector a, PolyVector b, PolyVector c){
	Polything arg[16];
	arg[0] = b.pos[0] - a.pos[0];
	arg[1] = b.pos[1] - a.pos[1];
	arg[2] = b.pos[2] - a.pos[2];
	arg[3] = b.pos[3] - a.pos[3];
	arg[4] = c.pos[0] - a.pos[0];
	arg[5] = c.pos[1] - a.pos[1];
	arg[6] = c.pos[2] - a.pos[2];
	arg[7] = c.pos[3] + a.pos[3];

	arg[8] = arg[0]*arg[0] +arg[1]*arg[1] +arg[2]*arg[2]; // ||b-a||2
	arg[9] = arg[0]*arg[4] +arg[1]*arg[5] +arg[2]*arg[6]; // c-a * b-a

	arg[10] =  arg[5]*arg[2]- arg[1]* arg[6];
	arg[11] =  arg[6]*arg[0]- arg[2]* arg[4];
	arg[12] =  arg[4]*arg[1]- arg[0]* arg[5];
	arg[13] =  arg[10]* arg[10] + arg[11] * arg[11] + arg[12] * arg[12];
	arg[14] = arg[3]* arg[3] - arg[8];
	arg[15] = arg[9] * arg[3] + arg[8] * arg[7];


	Polything f = arg[15]*arg[15] + arg[14]*arg[13];

	return(f.solve());
}


static double systemSolveSph(PolyVector a, PolyVector b, PolyVector c, PolyVector d){
	Polything arg[14];
	int i;
	for(i=0;i<3;i++){
	arg[i] = b.pos[i] - a.pos[i];
	arg[i+4] = c.pos[i] - a.pos[i];
	arg[i+8] = d.pos[i] - a.pos[i];
	}
	arg[3] = a.pos[3] - b.pos[3];
	arg[7] = a.pos[3] - c.pos[3];
	arg[11] = a.pos[3] + d.pos[3];

	arg[12] = arg[5] * arg[0] - arg[1] * arg[4];
	arg[13] = arg[6] * arg[0] - arg[2] * arg[4];
	arg[0] = arg[7] * arg[0] - arg[3] * arg[4]; //
	arg[4] = arg[6] * arg[1] - arg[2] * arg[5];
	arg[1] = arg[7] * arg[1] - arg[3] * arg[5]; //
	arg[5] = arg[7] * arg[2] - arg[3] * arg[6]; //


	arg[2] = arg[11] * arg[4] - arg[9] * arg[5] + arg[10] * arg[1]; // nx
	arg[6] = arg[8] * arg[5] - arg[11] * arg[13] + arg[10] * arg[0]; // ny
	arg[5] = arg[8] * arg[1] - arg[9] * arg[0] + arg[11] * arg[12]; // nx
	arg[1] = arg[8] * arg[4] - arg[9] * arg[13] + arg[10] * arg[12]; // nx

	arg[2].square();
	arg[6].square();
	arg[5].square();
	arg[1].square();

	Polything f = arg[2] + arg[6] + arg[5] - arg[1];
	return(f.solve());
}

static double systemSolveCyl(PolyVector a, PolyVector b, PolyVector c, PolyVector d){
	Polything arg[14];
	int i;
	for(i=0;i<3;i++){
	arg[i] = b.pos[i] - a.pos[i];
	arg[i+4] = d.pos[i] - c.pos[i];
	arg[i+8] = c.pos[i] - a.pos[i];
	}
	arg[3] = a.pos[3] - b.pos[3];
	arg[7] = c.pos[3] - d.pos[3];
	arg[11] = a.pos[3] + c.pos[3];

	arg[12] = arg[5] * arg[0] - arg[1] * arg[4];
	arg[13] = arg[6] * arg[0] - arg[2] * arg[4];
	arg[0] = arg[7] * arg[0] - arg[3] * arg[4]; //
	arg[4] = arg[6] * arg[1] - arg[2] * arg[5];
	arg[1] = arg[7] * arg[1] - arg[3] * arg[5]; //
	arg[5] = arg[7] * arg[2] - arg[3] * arg[6]; //


	arg[2] = arg[11] * arg[4] - arg[9] * arg[5] + arg[10] * arg[1]; // nx
	arg[6] = arg[8] * arg[5] - arg[11] * arg[13] + arg[10] * arg[0]; // ny
	arg[5] = arg[8] * arg[1] - arg[9] * arg[0] + arg[11] * arg[12]; // nx
	arg[1] = arg[8] * arg[4] - arg[9] * arg[13] + arg[10] * arg[12]; // nx

	arg[2].square();
	arg[6].square();
	arg[5].square();
	arg[1].square();

	Polything f = arg[2] + arg[6] + arg[5] - arg[1];
	return(f.solve());
}*/

/*
	int variint::eval(int time){
		int t;

		if (tpiv == bpiv) return(tpiv);
		else if (tpiv - bpiv > 0){
				if (tpiv < time) return(0x8FFFFFFF);
				*((uint32_t*)(&t)) = ((uint32_t)(time - bpiv)) * ((~((uint32_t) 0)) / ((uint32_t)(tpiv - bpiv)));
				return(t);
		}else{
				if (bpiv < time) return(0xFFFFFFFE);
				*((uint32_t*)(&t)) = ((uint32_t)(time - tpiv)) * ((~((uint32_t) 0)) / ((uint32_t)(bpiv - tpiv)));
				return(t);
		}
		}
*/
	int variint::eval(int time){
		int t;
		if (var == 0 ) return(piv);
		else {
				t = time - piv;
			if ((t >> 10) + 1 >= (0x7FFFFFFF / abs(var))) return( var > 0 ? 0x7FFFFFFF : 0x80000001);
		//	printf("%i\t%i\n", (t >> 10) * var, ((var * (t & 1023))>> 10) );
			if (t > 0) return( (t >> 10) * var + ((var * (t & 1023))>> 10) );
			else return( (t >> 10) * var + ((var * (1023-(t & 1023)))>> 10) );
		}
		}
	void variint::setValue(int value, int derivative, int time){
		var = derivative;
		if (var == 0) piv = value;
		else{
		piv = time - ((value / derivative)<< 10);
		}
		}





	/*
Polything::Polything(){memset(coef,'\0',sizeof(double)*(POLYTHING_DESIRED_ORDER+1));}
Polything::Polything(double *momments, double multiple, double exponent){
	int i = POLYTHING_DESIRED_ORDER-1;
	for(;i>=0;i--) coef[i] = momments[i];
	coef[POLYTHING_DESIRED_ORDER] = exp(exponent*timestep);
//	if ((exponent*timestep) > 2.0f){ // compute exact term (approximation is bad)
	for(i=1;i<=POLYTHING_DESIRED_ORDER;i++) coef[POLYTHING_DESIRED_ORDER] = (coef[POLYTHING_DESIRED_ORDER] - 1) * i / (exponent*timestep);
//	}
	coef[POLYTHING_DESIRED_ORDER] *= (multiple * pow(exponent,(double)POLYTHING_DESIRED_ORDER) / POLYTHING_DESIRED_ORDER_FACTORIAL);
}


void Polything::clear(){memset(coef,'\0',sizeof(double)*(POLYTHING_DESIRED_ORDER+1));}

void Polything::ODEextends(ODEfunctions ode, double factor, double exp){
	int i,j;
	double tmp;
	switch(ode){
		case LFHP_ODE_2ND_COS: // the 2 first momments are used
				for(i = 2;i<POLYTHING_DESIRED_ORDER;i++){
					tmp =coef[0];
					for(j=1;j<=i;j++){
						if (((i-j) % 2) == 1) tmp /= j;
						else tmp = (tmp + coef[j] * ((((i-j) % 4)  == 0)? 1 : -1)  ) /j;
					}
					coef[i] = factor * tmp / ((i-1)*(i));
				}
			break;
		case LFHP_ODE_2ND_SIN: // the 2 first momments are used
				// motion on a cylinder, in angular coordonates
				for(i = 2;i<POLYTHING_DESIRED_ORDER;i++){
					tmp =0;
					for(j=1;j<=i;j++){
						if (((i-j) % 2) == 0) tmp /= j;
						else tmp = (tmp + coef[j] * ((((i-j) % 4)  == 1)? 1 : -1)  ) /j;
					}
					coef[i] = factor * tmp / ((i-1)*(i));
				}
			break;
	}
}

Polything& Polything::operator+=(Polything const & other) {
	int i = POLYTHING_DESIRED_ORDER-1;
	for(;i>=0;i--) coef[i] += other.coef[i];
	coef[POLYTHING_DESIRED_ORDER] += other.coef[POLYTHING_DESIRED_ORDER];
	return(*this);
}

Polything& Polything::operator-=(Polything const & other) {
	int i = POLYTHING_DESIRED_ORDER-1;
	for(;i>=0;i--) coef[i] -= other.coef[i];
	coef[POLYTHING_DESIRED_ORDER] += other.coef[POLYTHING_DESIRED_ORDER];
	return(*this);
}

Polything& Polything::operator*=(Polything const & other) {
	int i,j;
	double sum = 0.0f;
	double multi =0.0f;

	for(i= 0;i<POLYTHING_DESIRED_ORDER+1;i++) {
		multi += fabs(coef[POLYTHING_DESIRED_ORDER-i]);
		sum += fabs(other.coef[i]) * multi;
		multi *= timestep;
	}
	coef[POLYTHING_DESIRED_ORDER] = sum;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		sum =0;
		for(j= 0;j<=i;j++) {
			sum += coef[i-j] * other.coef[j];
		}
		coef[i] = sum;
	}
	return(*this);
}*/
/*
Polything& Polything::operator*(Polything const & other, int const & fact) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] = other.coef[i] * fact;
	}
	return(*this);
}

Polything& Polything::operator*(Polything const & other, double const & fact) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] = other.coef[i] * fact;
	}
	return(*this);
}

Polything& Polything::operator/(Polything const & other, int const & fact) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] = other.coef[i] / fact;
	}
	return(*this);
}

Polything& Polything::operator/(Polything const & other, double const & fact) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] = other.coef[i] / fact;
	}
	return(*this);
}
*/
/*
Polything Polything::operator*(int const & other) const{
	int i;
	Polything newpoly;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		newpoly.coef[i] = coef[i] * other;
	}
	return(newpoly);
}

Polything Polything::operator*(double const & other) const{
	int i;
	Polything newpoly;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		newpoly.coef[i] = coef[i] * other;
	}
	return(newpoly);
}

Polything Polything::operator/(int const & other) const{
	int i;
	Polything newpoly;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		newpoly.coef[i] = coef[i] / other;
	}
	return(newpoly);
}

Polything Polything::operator/(double const & other) const{
	int i;
	Polything newpoly;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		newpoly.coef[i] = coef[i] / other;
	}
	return(newpoly);
}

Polything& Polything::operator+=(int const & other) {coef[0] += other;return(*this);}
Polything& Polything::operator+=(double const & other) {coef[0] += other;return(*this);}
Polything& Polything::operator-=(int const & other) {coef[0] -= other;return(*this);}
Polything& Polything::operator-=(double const & other) {coef[0] -= other;return(*this);}
Polything Polything::operator+(int const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	newpoly.coef[0] = coef[0] + other;
	return(newpoly);
}
Polything Polything::operator+(double const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	newpoly.coef[0] = coef[0] + other;
	return(newpoly);
}
Polything Polything::operator-(int const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	newpoly.coef[0] = coef[0] - other;
	return(newpoly);
}
Polything Polything::operator-(double const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	newpoly.coef[0] = coef[0] - other;
	return(newpoly);
}

Polything& Polything::operator*=(int const & other) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] *= other;
	}
	return(*this);
}

Polything& Polything::operator*=(double const & other) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] *= other;
	}
	return(*this);
}

Polything& Polything::operator/=(int const & other) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] /= other;
	}
	return(*this);
}

Polything& Polything::operator/=(double const & other) {
	int i;
	for(i= POLYTHING_DESIRED_ORDER-1;i>=0;i--) {
		coef[i] /= other;
	}
	return(*this);
}

Polything Polything::operator+(Polything const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	int i;
	for(i= POLYTHING_DESIRED_ORDER;i>=0;i--) {
	newpoly.coef[i] = coef[i] + other.coef[i];
	}
	return(newpoly);
}
Polything Polything::operator-(Polything const & other) const{
	Polything newpoly;
	memcpy(newpoly.coef,coef,sizeof(double)* (POLYTHING_DESIRED_ORDER+1));
	int i;
	for(i= 0;i<POLYTHING_DESIRED_ORDER;i++) {
	newpoly.coef[i] = coef[i] - other.coef[i];
	}
	newpoly.coef[i] = coef[i] + other.coef[i];
	return(newpoly);
}
Polything Polything::operator*(Polything const & other) const{
	int i,j;
	Polything newpoly;
	for(i= 0;i<POLYTHING_DESIRED_ORDER;i++) {
		newpoly.coef[i] = coef[i] * other.coef[0];
		for(j=1 ;j<=i;j++){
		newpoly.coef[i] += coef[i-j] * other.coef[j];
		}
	}
	newpoly.coef[i] = fabs(coef[i] * other.coef[0]);
	for(j=1 ;j<=i;j++){
		newpoly.coef[i] += fabs(coef[i-j] * other.coef[j]);
	}
	return(newpoly);
}

Polything Polything::getsquare(){
	int i,j;
	Polything newpoly;
	for(i= 0;i<POLYTHING_DESIRED_ORDER;i++) {
		newpoly.coef[i] = coef[i] * coef[0];
		for(j=1 ;j<=i;j++){
		newpoly.coef[i] += coef[i-j] * coef[j];
		}
	}
	newpoly.coef[i] = fabs(coef[i] * coef[0]);
	for(j=1 ;j<=i;j++){
		newpoly.coef[i] += fabs(coef[i-j] * coef[j]);
	}
	return(newpoly);
}

void Polything::square(){
	int i,j;
	double tmp[POLYTHING_DESIRED_ORDER+1];
	for(i= 0;i<POLYTHING_DESIRED_ORDER;i++) {
		tmp[i] = coef[i] * coef[0];
		for(j=1 ;j<=i;j++){
		tmp[i] += coef[i-j] * coef[j];
		}
	}
	tmp[i] = fabs(coef[i] * coef[0]);
	for(j=1 ;j<=i;j++){
		tmp[i] += fabs(coef[i-j] * coef[j]);
	}
	memcpy(coef,tmp,sizeof(double)*(POLYTHING_DESIRED_ORDER+1));
}



double Polything::quadsolve(){return(quadsolve(coef));}
double Polything::orientedquadsolve(){return(orientedquadsolve(coef));}
double Polything::negorientedquadsolve(){return(negorientedquadsolve(coef));}


double Polything::solve(){
	double out = timestep;
	double tmp = coef[0];
	double tmp2;
	double tmp3;
	double pivot;
	if (coef[0] <0) {return(out);}
	int i,j;
	for(i=1;i<=POLYTHING_DESIRED_ORDER;i++){
		tmp = tmp / out + coef[i];
		if (tmp < 0) break;
	}
	if (i > POLYTHING_DESIRED_ORDER) return(out); // no zero
	pivot = out * (1.0f - tmp / coef[i]);

	if (coef[1] >= 0){
	// high order approach
		if (coef[3] >= 0) {
			tmp = coef[2] + coef[1] / pivot;
			tmp2 = coef[4] + coef[3] / pivot;
		}else{
			tmp = tmp + coef[1] / pivot + coef[3] * pivot;
			tmp2 = coef[4];
		}
			tmp3 = tmp * tmp - 4 * coef[0] * tmp2;
			if (tmp3< 0) return(pivot);
			tmp3 = sqrt(tmp3);
			tmp3 = (tmp2 > 0)^(tmp3 < fabs(tmp)) ? -tmp - tmp3 : -tmp + tmp3;
			return(sqrt( tmp3 / (2*tmp2)));

	}else{
	// low order approach
	tmp = coef[POLYTHING_DESIRED_ORDER];
	tmp2 =  coef[POLYTHING_DESIRED_ORDER-1];
	for(i = POLYTHING_DESIRED_ORDER-1;i>=3;i--){
		if (tmp < 0) {
		tmp = tmp * pivot + tmp2;
		tmp2 = coef[i-1];
		}else{
		tmp2 = coef[i-1] + pivot*pivot*tmp;
		tmp = tmp2 - 2 * pivot*tmp;
		}
	}
	tmp3 = tmp * tmp - 4 * coef[0] * tmp2;
	if (tmp3< 0) return(pivot);
	tmp3 = sqrt(tmp3);
	tmp3 = (tmp2 > 0)^(tmp3 < fabs(tmp)) ? -tmp - tmp3 : -tmp + tmp3;
	return(tmp3);
	}


}
*/

/*

double Polything::solve(){
	double out = timestep;


}
*/

/*
//	pivot = tmp;
	j =i;
	min = tmp;
//	pivot = timestep*coef[0] /(coef[0] - min * pow(timestep, POLYTHING_DESIRED_ORDER -i)  );
	for(i++;i<=POLYTHING_DESIRED_ORDER;i++){
		tmp = tmp / out + coef[i]
		if (tmp < min/out) min = tmp;
		else min = min/out;
	}
	min = min * pow(timestep,j);

	pivot = timestep*coef[0] /(coef[0] - pivot); // new pivot point


*/
/*double Polything::solve(){
	int ite=POLYTHING_DESIRED_ITERATION_NB;
	double tmp[3];
	tmp[0] = coef[0];
	double zero = timestep;
	int i;
	bool positive = (tmp[0]>0);
	bool extends= true;
	for(i =1;i<POLYTHING_DESIRED_ORDER+1;i++){
		if (i == POLYTHING_DESIRED_ORDER){
			tmp[1] = fabs(tmp[0]) / coef[i];
			if ((tmp[1] > zero)&&(( ( ((tmp[1] - zero) / tmp[1]) < 0.02f) ) ||(ite == 1))  ) return(zero);
			//		zero = (zero + fabs(tmp) / coef[i]) / (2); // coef[i] cannot be zero here
			//		zero = (zero*ite + fabs(tmp) / coef[i]) / (1+ite); // coef[i] cannot be zero here
			if ((tmp[1] > zero/2)&&(tmp[1] / 2 < zero)) zero = tmp[1]; // coef[i] cannot be zero here
			else if (tmp[1] > zero) zero =zero*2;
			else zero =zero/2;
			i =0;
			tmp[0] = coef[0];
			if (ite >1) ite--;
		}else{
			if ((((tmp[1] =  coef[i] + tmp[0]/zero) < 0)^(positive))&&(tmp[1] != 0.0f)) {
				tmp[0] = tmp[1];
			}else if ((i+1 < POLYTHING_DESIRED_ORDER)&&(((tmp[1] = coef[i+1] - coef[i]*coef[i]/ (4 * tmp[0])) < 0)^(positive))&&(tmp[1] != 0.0f)){
				tmp[0] = tmp[1]; i=i+1;
			}else break;
		}
	}

	if (i == POLYTHING_DESIRED_ORDER+1) return(zero);
	zero = pow(-coef[0] / coef[i],1.0f/i); // approximate to the closest zero
	i= POLYTHING_DESIRED_ORDER;
	while (coef[i] == 0) i--;
	tmp[2] = coef[i];
	tmp[1] = coef[i-1];
	for(;i>2;i--){
		if ((coef[i] < 0)^(positive)){
			tmp[2] = tmp[1] + tmp[2] * zero;
			tmp[1] = coef[i-2];
		}else{
			tmp[0] = tmp[2];
			tmp[2] = tmp[1] + tmp[2] * zero * 2;
			tmp[1] = coef[i-2] - zero * zero * tmp[0];
			extends = false;
		}
	}
	tmp[0] = coef[0];
	tmp[0] = quadsolve(tmp);
	if ((tmp[0]<zero)||(extends)) return(tmp[0]);
	else return(zero);
}*/
/*
double Polything::eval(double const time) const{
	double out = coef[POLYTHING_DESIRED_ORDER];
	int i = POLYTHING_DESIRED_ORDER-1;
	for(;i>=0;i--){
		out = coef[i] + out * time;
	}
	return(out);
}

double Polything::evalDerivative(double const time) const{
	int i = POLYTHING_DESIRED_ORDER-1;
	double out = coef[POLYTHING_DESIRED_ORDER] * i;
	for(i--;i>=1;i--){
		out = coef[i] + out * time * i;
	}
	return(out);
}

void Polything::show(FILE* f){
	fprintf(f,"Polynominal: ");
	int i;
	int step = 1000;
	fprintf(f,"%.2f %c %.2fx ",coef[0], (coef[1]>=0) ? '+' : '-', fabs(coef[1]));
	for(i=2;i<POLYTHING_DESIRED_ORDER;i++) fprintf(f,"%c %.2f x%i ", (coef[i]>=0) ? '+' : '-' ,fabs(coef[i]),i);
	fprintf(f,"+- %fx%i ",coef[POLYTHING_DESIRED_ORDER],i);
	fprintf(f," (zero = %f with timestep %.2f)\n",solve(),timestep);
 }

Polything PolyVector::dotProduct(double const * const & other) const{
Polything newpoly;
int i;
memset(newpoly.coef,'\0',sizeof(double) *(POLYTHING_DESIRED_ORDER+1));
for(i=0;i<order;i++){
	newpoly.coef[i] = pos[0].coef[i] * other[0]
		            + pos[1].coef[i] * other[1]
					+ pos[2].coef[i] * other[2];
}
return(newpoly);
}

PolyVector PolyVector::crossProduct(double const * const & other) const{
PolyVector newpoly;
newpoly.order = order;
int i;
for(i=0;i<order;i++){
	newpoly.pos[0].coef[i] = pos[1].coef[i] * other[2] - pos[2].coef[i] * other[1];
	newpoly.pos[1].coef[i] = pos[2].coef[i] * other[0] - pos[0].coef[i] * other[2];
	newpoly.pos[2].coef[i] = pos[0].coef[i] * other[1] - pos[1].coef[i] * other[0];
}
return(newpoly);
}

Polything PolyVector::squarrenorm() const{
Polything newpoly;
int i,j;
memset(newpoly.coef,'\0',sizeof(double) *(POLYTHING_DESIRED_ORDER+1));
for(i=0;i<order;i++){
	for(j=i;(j<order)&&(j+i<POLYTHING_DESIRED_ORDER);j++){
		if (i == j) {
			newpoly.coef[2*i] += pos[0].coef[i] *pos[0].coef[i] +pos[1].coef[i] *pos[1].coef[i] +pos[2].coef[i] *pos[2].coef[i];
		} else{
			newpoly.coef[i+j] += 2*(pos[0].coef[i] *pos[0].coef[j] +pos[1].coef[i] *pos[1].coef[j] +pos[2].coef[i] *pos[2].coef[j]);
		}
	}
}
return(newpoly);
}


void PolyVector::eval(double const time, double * const out) const{
int i = order-1;
out[0] = pos[0].coef[i];
out[1] = pos[1].coef[i];
out[2] = pos[2].coef[i];
for(i--;i>=0;i--){
out[0] = pos[0].coef[i] + time * out[0];
out[1] = pos[1].coef[i] + time * out[1];
out[2] = pos[2].coef[i] + time * out[2];
}
if (_isnan(*out)){
int asdjfkas=34829;
}
}

void PolyVector::evalDerivative(double const time, double * const out) const{
int i = order-1;
out[0] = pos[0].coef[i]*i;
out[1] = pos[1].coef[i]*i;
out[2] = pos[2].coef[i]*i;
for(i--;i>=1;i--){
out[0] = pos[0].coef[i] + time * out[0]*i;
out[1] = pos[1].coef[i] + time * out[1]*i;
out[2] = pos[2].coef[i] + time * out[2]*i;
}
}

int PolyVector::operator==(PolyVector const & other) const{
	bool max = order > other.order;
	int k = (max) ? other.order : order;
	int i,j;
	for(i=0;i<k;i++)
		for(j=0;j<4;j++)
			if (pos[j].coef[i] != other.pos[j].coef[i]) return(false);
	if (max){
		for(;i<order;i++)
			for(j=0;j<4;j++)
				if (pos[j].coef[i] != 0.0f) return(false);
	}else{
		for(;i<other.order;i++)
			for(j=0;j<4;j++)
				if (other.pos[j].coef[i] != 0.0f) return(false);
	}
	return(true);
}

PolyVector PolyVector::operator-(double const * const & other) const{
PolyVector newpoly;
int i;
	newpoly.order = order;
	newpoly.pos[0].coef[0] = pos[0].coef[0] - other[0];
	newpoly.pos[1].coef[0] = pos[1].coef[0] - other[1];
	newpoly.pos[2].coef[0] = pos[2].coef[0] - other[2];
	newpoly.pos[3].coef[0] = pos[3].coef[0] - other[3];

for(i=1;i<order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i];
}

return(newpoly);
}

PolyVector PolyVector::operator+(PolyVector const & other) const{
	PolyVector newpoly;
	int i;
	if (order == other.order) {
	newpoly.pos[0] = pos[0] + other.pos[0];
	newpoly.pos[1] = pos[1] + other.pos[1];
	newpoly.pos[2] = pos[2] + other.pos[2];
	newpoly.pos[3] = pos[3] + other.pos[3];
	newpoly.order = order;
	}else if (order > other.order){
	newpoly.order = order;
	for(i=0;i<other.order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i] + other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i] + other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i] + other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i] + other.pos[3].coef[i];
	}
	for(;i<order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i];
	}
	}else{
	newpoly.order = other.order;
	for(i=0;i<order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i] + other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i] + other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i] + other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i] + other.pos[3].coef[i];
	}
	for(;i<other.order;i++){
	newpoly.pos[0].coef[i] = other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = other.pos[3].coef[i];
	}
	}
	return(newpoly);
}


PolyVector PolyVector::operator-(PolyVector const & other) const{
	PolyVector newpoly;
	int i;
	if (order == other.order) {
	newpoly.pos[0] = pos[0] - other.pos[0];
	newpoly.pos[1] = pos[1] - other.pos[1];
	newpoly.pos[2] = pos[2] - other.pos[2];
	newpoly.pos[3] = pos[3] - other.pos[3];
	newpoly.order = order;
	}else if (order > other.order){
	newpoly.order = order;
	for(i=0;i<other.order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i] - other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i] - other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i] - other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i] - other.pos[3].coef[i];
	}
	for(;i<order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i];
	}
	}else{
	newpoly.order = other.order;
	for(i=0;i<order;i++){
	newpoly.pos[0].coef[i] = pos[0].coef[i] - other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = pos[1].coef[i] - other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = pos[2].coef[i] - other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = pos[3].coef[i] - other.pos[3].coef[i];
	}
	for(;i<other.order;i++){
	newpoly.pos[0].coef[i] = -other.pos[0].coef[i];
	newpoly.pos[1].coef[i] = -other.pos[1].coef[i];
	newpoly.pos[2].coef[i] = -other.pos[2].coef[i];
	newpoly.pos[3].coef[i] = -other.pos[3].coef[i];
	}
	}
	return(newpoly);
}*/
/*
PolyVector& PolyVector::operator-(PolyVector const & other) const{
	PolyVector newpoly;
	newpoly.pos[0] = pos[0] - other.pos[0];
	newpoly.pos[1] = pos[1] - other.pos[1];
	newpoly.pos[2] = pos[2] - other.pos[2];
	newpoly.pos[3] = pos[3] - other.pos[3];
	newpoly.order = order;
	return(newpoly);
}*/

/*
PolyVector PolyVector::operator/(double const & other) const{
	PolyVector newpoly;
	newpoly.pos[0] = pos[0] / other;
	newpoly.pos[1] = pos[1] / other;
	newpoly.pos[2] = pos[2] / other;
	newpoly.pos[3] = pos[3] / other;
	newpoly.order = order;
	return(newpoly);
}

// EZ Matrix
PolyMatrix::PolyMatrix(PolyVector &centeroffset, SQuaternion ori){
	int tmpint[10];
	tmpint[0] = (int)ori.x;
	tmpint[1] = (int)ori.y;
	tmpint[2] = (int)ori.z;
	tmpint[3] = -tmpint[0] * ori.w;
	tmpint[4] = -tmpint[1] * ori.w;
	tmpint[5] = -tmpint[2] * ori.w;
	tmpint[6] = tmpint[1] * tmpint[2];
	tmpint[7] = tmpint[0] * tmpint[2];
	tmpint[8] = tmpint[0] * tmpint[1];

	tmpint[9] = tmpint[2] * tmpint[2];
	tmpint[2] = tmpint[1] * tmpint[1];
	tmpint[1] = tmpint[9];
	tmpint[9] = tmpint[0] * tmpint[0];
	tmpint[0] = tmpint[1] + tmpint[2];
	tmpint[1] = tmpint[9] + tmpint[1];
	tmpint[2] = tmpint[9] + tmpint[2];

	memset(pos+3, '\0', sizeof(Polything)*9);
	pos[3].coef[0] = (536870912 - tmpint[0]) / 536870912.0f;
	pos[4].coef[0] = (tmpint[8] +tmpint[5]) / 536870912.0f;
	pos[5].coef[0] = (tmpint[7] -tmpint[4]) / 536870912.0f;
	pos[6].coef[0] = (tmpint[8] -tmpint[5]) / 536870912.0f;
	pos[7].coef[0] = (536870912 - tmpint[1]) / 536870912.0f;
	pos[8].coef[0] = (tmpint[6] -tmpint[3]) / 536870912.0f;
	pos[9].coef[0] = (tmpint[7] +tmpint[4]) / 536870912.0f;
	pos[10].coef[0] = (tmpint[6] +tmpint[3]) / 536870912.0f;
	pos[11].coef[0] = (536870912 - tmpint[2]) / 536870912.0f;

	pos[0] = centeroffset.pos[0];
	pos[1] = centeroffset.pos[1];
	pos[2] = centeroffset.pos[2];
	order = centeroffset.order;
}

// Da purely EVIL polyMatrix!
PolyMatrix::PolyMatrix(PolyVector &rotapos, PolyVector &centeroffset, SQuaternion ori, PolyVector &angular){

	// angular should have coef[0] = 0 for x,y and z
	Polything N = angular.squarrenorm();
	Polything CS = N.getsquare();
	Polything SS = (N + (CS / 5))* (2.0f/(3.0f*268435456.0f)); // (1 - 2/3 N  + 2/15 N^2)
	CS = (N + CS * (2.0f /15.0f))/(3.0f *268435456.0f);	// (1 - 1/3 N  + 2/45 N^2)
	CS.coef[0] = 1.0f / 268435456.0f;
	SS.coef[0] = 1.0f / 268435456.0f;

	Polything tmp[3];
	tmp[0] = angular.pos[0].getsquare();
	tmp[1] = angular.pos[1].getsquare();
	tmp[2] = angular.pos[2].getsquare();

	int tmpint[10];
	tmpint[0] = (int)ori.x;
	tmpint[1] = (int)ori.y;
	tmpint[2] = (int)ori.z;
	tmpint[3] = -tmpint[0] * ori.w;
	tmpint[4] = -tmpint[1] * ori.w;
	tmpint[5] = -tmpint[2] * ori.w;
	tmpint[6] = tmpint[1] * tmpint[2];
	tmpint[7] = tmpint[0] * tmpint[2];
	tmpint[8] = tmpint[0] * tmpint[1];

	tmpint[9] = tmpint[2] * tmpint[2];
	tmpint[2] = tmpint[1] * tmpint[1];
	tmpint[1] = tmpint[9];
	tmpint[9] = tmpint[0] * tmpint[0];
	tmpint[0] = tmpint[1] + tmpint[2];
	tmpint[1] = tmpint[9] + tmpint[1];
	tmpint[2] = tmpint[9] + tmpint[2];


*/

/*
	pos[3] = SS*(tmpint[0]*(tmp[1] + tmp[2])+ angular.pos[0]*
		         ((tmpint[2]-tmpint[5])*angular.pos[1]+ (tmpint[1]+tmpint[4])*angular.pos[2]))
		   + CS*(((tmpint[2]-tmpint[5])*angular.pos[2]+ (tmpint[1]+tmpint[4])*angular.pos[1]));

	pos[4] =



*/

/*
	pos[3] = SS*((tmp[1] + tmp[2])* -1.0f);

	pos[7] = SS*((tmp[0] + tmp[2])* -1.0f);
	pos[7].coef[0] += 1.0f;
	pos[11] = SS*((tmp[0] + tmp[1])* -1.0f);
	pos[11].coef[0] += 1.0f;

	tmp[0] = angular.pos[0]* angular.pos[1]*SS;
	tmp[1] = angular.pos[2]*CS;
	pos[4] = tmp[0] + tmp[1];
	pos[6] = tmp[0] - tmp[1];

	tmp[0] = angular.pos[0]* angular.pos[2]*SS;
	tmp[1] = angular.pos[1]*CS;
	pos[5] = tmp[0] - tmp[1];
	pos[9] = tmp[0] + tmp[1];

	tmp[0] = angular.pos[1]* angular.pos[2]*SS;
	tmp[1] = angular.pos[0]*CS;
	pos[8] = tmp[0] + tmp[1];
	pos[10]= tmp[0] - tmp[1];

	pos[3].coef[0] = (536870912 - tmpint[0]) / 536870912.0f;
	pos[4].coef[0] = (tmpint[8] +tmpint[5]) / 536870912.0f;
	pos[5].coef[0] = (tmpint[7] -tmpint[4]) / 536870912.0f;
	pos[6].coef[0] = (tmpint[8] -tmpint[5]) / 536870912.0f;
	pos[7].coef[0] = (536870912 - tmpint[1]) / 536870912.0f;
	pos[8].coef[0] = (tmpint[6] -tmpint[3]) / 536870912.0f;
	pos[9].coef[0] = (tmpint[7] +tmpint[4]) / 536870912.0f;
	pos[10].coef[0] = (tmpint[6] +tmpint[3]) / 536870912.0f;
	pos[11].coef[0] = (536870912 - tmpint[2]) / 536870912.0f;
}

PolyVector PolyMatrix::operator*(PolyVector const & traj) const{
	PolyVector newpoly;
	int i,j;
	PolyVector tmp = traj;
	tmp.pos[0] -= pos[0];
	tmp.pos[1] -= pos[1];
	tmp.pos[2] -= pos[2];
//	if (traj.order == order){
//	if (pos[3].coef[0] != 1.0f){
//	int wanasee=1;
//	}
	newpoly.pos[0] = pos[3] * tmp.pos[0] + pos[6] * tmp.pos[1] + pos[9] * tmp.pos[2];
	newpoly.pos[1] = pos[4] * tmp.pos[0] + pos[7] * tmp.pos[1] + pos[10] * tmp.pos[2];
	newpoly.pos[2] = pos[5] * tmp.pos[0] + pos[8] * tmp.pos[1] + pos[11] * tmp.pos[2];
	newpoly.pos[3] = traj.pos[3];
	newpoly.order =  traj.order > order ? traj.order :order;
*/

/*
	}else if (order > traj.order){

	PolyVector tmp = traj;
	for(i = traj.order;i<order;i++){
	for(j=0;j<4;j++) tmp.pos[j].coef[i] =0.0f;
	}

	newpoly.pos[0] = pos[3] * tmp.pos[0] + pos[6] * tmp.pos[1] + pos[9] * tmp.pos[2] - pos[0];
	newpoly.pos[1] = pos[4] * tmp.pos[0] + pos[7] * tmp.pos[1] + pos[10] * tmp.pos[2] - pos[1];
	newpoly.pos[2] = pos[5] * tmp.pos[0] + pos[8] * tmp.pos[1] + pos[11] * tmp.pos[2] - pos[2];
	newpoly.pos[3] = tmp.pos[3];
	newpoly.order = order;

	}else{


	}
	return(newpoly);
}*/

/*
Event::Event(): RBTreeNode<int,IntComparator>(){}
Event::~Event(){}

int Event::getIndex(){return(getTime());}

bool Event::autoDelete(){return(true);}

PriorityQueue::PriorityQueue(): sema(false){
	time.speed = 0.0f;
	time.offset = 0.0f;
}


void PriorityQueue::Insert(Event* ev){
toinsert.push_back(ev);
}

void PriorityQueue::start(){
	time.speed = 1.0f;
	time.offset = time.offset + clock() * time.speed;
}
void PriorityQueue::pause(){
	time.offset = time.offset - clock() * time.speed;
	time.speed = 0.0f;
}

void PriorityQueue::dotilldone(){
	int toptime;
	int i;
	time();
	int tmptime = time.time;
//	while(1){
//	while( (first != NULL) && ( (time.time > (toptime = first->index)) || ( toptime < (time.time = (int)(clock() * speed - offset))) ) ) {

	if (( (i = toinsert.size()) > 0)&&(!sema)){
		sema = true;
		for(i--;i>=0;i--){list.Insert(toinsert[i]);}
		toinsert.clear();
		sema = false;
	}

	while( (list.first != NULL) && (tmptime > (toptime = list.first->getIndex())) ) {
			time.time = toptime;
	//		fprintf(Dbg,"poping:\n");
			Event* hoho = (Event*) list.Pop();
			(*hoho)();
			if (hoho->autoDelete()) delete(hoho);
			if (( (i = toinsert.size()) > 0)&&(!sema)){
				sema = true;
				for(i--;i>=0;i--){list.Insert(toinsert[i]);}
				toinsert.clear();
				sema = false;
			}
	}

//	time();
//	tmptime = time.time;
//	if (toptime > tmptime) break;
//	}
}

void PriorityQueue::do_one(){
	int i;
	if (( (i = toinsert.size()) > 0)&&(!sema)){
		sema = true;
		for(i--;i>=0;i--){list.Insert(toinsert[i]);}
		toinsert.clear();
		sema = false;
	}

	if (list.first != NULL){
	Event* hoho = (Event*) list.Pop();
		(*hoho)();
		if (hoho->autoDelete()) delete(hoho);
	}
	}
 code generator
		int i;
	int i;
	for(i = 0; i<16;i++){
		printf("case '\\x%x': return(value >> %i); case '\\x%x': return(-(value >> %i)); case '\\x%x': return(-((-value) >> %i)); case '\\x%x': return((-value) >> %i);\n", i * 8,15-i, i * 8 + 128,15-i, i * 8 + 256,15-i, i * 8 + 384,15-i);
		printf("case '\\x%x': return((value >> %i) + (value >> %i)); case '\\x%x': return(-(value >> %i) - (value >> %i)); case '\\x%x': return(-((-value) >> %i) - ((-value) >> %i)); case '\\x%x': return(((-value) >> %i) + ((-value) >> %i));\n", i * 8+1,15-i,18-i, i * 8 + 129,15-i,18-i, i * 8 + 257,15-i,18-i, i * 8 + 385,15-i,18-i);
		printf("case '\\x%x': return((value >> %i) + (value >> %i)); case '\\x%x': return(-(value >> %i) - (value >> %i)); case '\\x%x': return(-((-value) >> %i) - ((-value) >> %i)); case '\\x%x': return(((-value) >> %i) + ((-value) >> %i));\n", i * 8+2,15-i,17-i, i * 8 + 130,15-i,17-i, i * 8 + 258,15-i,17-i, i * 8 + 386,15-i,17-i);
		printf("case '\\x%x': return(((value >> %i) - (value >> %i)) - (value >> %i)); case '\\x%x': return(((value >> %i) - (value >> %i)) + (value >> %i) ); case '\\x%x': return(((value >> %i) - (value >> %i)) + ((-value) >> %i) ); case '\\x%x': return(((value >> %i) - (value >> %i)) - ((-value) >> %i));\n", i * 8+3,14-i,16-i,18-i, i * 8 + 131,16-i,14-i,18-i, i * 8 + 259,14-i,16-i,18-i, i * 8 + 387,16-i,14-i,18-i);
		printf("case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i));\n", i * 8+4,14-i,16-i, i * 8 + 132,16-i,14-i, i * 8 + 260,14-i,16-i, i * 8 + 388,16-i,14-i);
		printf("case '\\x%x': return(((value >> %i) - (value >> %i)) + (value >> %i)); case '\\x%x': return(((value >> %i) - (value >> %i)) + (value >> %i) ); case '\\x%x': return(((value >> %i) - (value >> %i)) + ((-value) >> %i) ); case '\\x%x': return(((value >> %i) - (value >> %i)) - ((-value) >> %i));\n", i * 8+5,14-i,16-i,18-i, i * 8 + 133,16-i,14-i,18-i, i * 8 + 261,14-i,16-i,18-i, i * 8 + 389,16-i,14-i,18-i);
		printf("case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i));\n", i * 8+6,14-i,17-i, i * 8 + 134,17-i,14-i, i * 8 + 262,14-i,17-i, i * 8 + 390,17-i,14-i);
		printf("case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i)); case '\\x%x': return((value >> %i) - (value >> %i));\n", i * 8+7,14-i,18-i, i * 8 + 135,18-i,14-i, i * 8 + 263,14-i,18-i, i * 8 + 391,18-i,14-i);
		}

*/



	/*

void StateAutomata_State::randomize(){

	int i;
	for(i=0;i<176;i++) nextarg_weights[i] = rand() & 0x000000FF;
	for(i=0;i<22;i++) weights[i] = rand() & 0x000000FF;
	selector = rand() ^ (rand() << 12) ^ (rand() << 24);
	for(i=0;i<4;i++) nextstate[i] = rand() ^ (rand() << 12) ^ (rand() << 24);
	for(i=0;i<8;i++) queries[i] = rand() ^ (rand() << 12) ^ (rand() << 24);
	}

int StateAutomata_State::app_weight(char weight, int value){

	switch((((uint32_t)weight) & 0x000000FF) | ((value < 0) ? 256 : 0) ){
case 0x0: return(0); case 0x80: return(0); case 0x100: return(0); case 0x180: return(0);
case 0x1: return(value >> 17); case 0x81: return(-(value >> 17)); case 0x101: return(-((-value) >> 17)); case 0x181: return((-value) >> 17);
case 0x2: return(value >> 16); case 0x82: return(-(value >> 16)); case 0x102: return(-((-value) >> 16)); case 0x182: return((-value) >> 16);
case 0x3: return((value >> 15) - (value >> 17)); case 0x83: return((value >> 17) - (value >> 15)); case 0x103: return((value >> 15) - (value >> 17)); case 0x183: return((value >> 17) - (value >> 15));
case 0x4: return(value >> 15); case 0x84: return(-(value >> 15)); case 0x104: return(-((-value) >> 15)); case 0x184: return((-value) >> 15);
case 0x5: return((value >> 15) + (value >> 17)); case 0x85: return(-(value >> 15) - (value >> 17)); case 0x105: return(-((-value) >> 15) - ((-value) >> 17)); case 0x185: return(((-value) >> 15) + ((-value) >> 17));
case 0x6: return((value >> 14) - (value >> 16)); case 0x86: return((value >> 16) - (value >> 14)); case 0x106: return((value >> 14) - (value >> 16)); case 0x186: return((value >> 16) - (value >> 14));
case 0x7: return((value >> 14) - (value >> 17)); case 0x87: return((value >> 17) - (value >> 14)); case 0x107: return((value >> 14) - (value >> 17)); case 0x187: return((value >> 17) - (value >> 14));
case 0x8: return(value >> 14); case 0x88: return(-(value >> 14)); case 0x108: return(-((-value) >> 14)); case 0x188: return((-value) >> 14);
case 0x9: return((value >> 14) + (value >> 17)); case 0x89: return(-(value >> 14) - (value >> 17)); case 0x109: return(-((-value) >> 14) - ((-value) >> 17)); case 0x189: return(((-value) >> 14) + ((-value) >> 17));
case 0xa: return((value >> 14) + (value >> 16)); case 0x8a: return(-(value >> 14) - (value >> 16)); case 0x10a: return(-((-value) >> 14) - ((-value) >> 16)); case 0x18a: return(((-value) >> 14) + ((-value) >> 16));
case 0xb: return(((value >> 13) - (value >> 15)) - (value >> 17)); case 0x8b: return(((value >> 15) - (value >> 13)) + (value >> 17) ); case 0x10b: return(((value >> 13) - (value >> 15)) + ((-value) >> 17) ); case 0x18b: return(((value >> 15) - (value >> 13)) - ((-value) >> 17));
case 0xc: return((value >> 13) - (value >> 15)); case 0x8c: return((value >> 15) - (value >> 13)); case 0x10c: return((value >> 13) - (value >> 15)); case 0x18c: return((value >> 15) - (value >> 13));
case 0xd: return(((value >> 13) - (value >> 15)) + (value >> 17)); case 0x8d: return(((value >> 15) - (value >> 13)) + (value >> 17) ); case 0x10d: return(((value >> 13) - (value >> 15)) + ((-value) >> 17) ); case 0x18d: return(((value >> 15) - (value >> 13)) - ((-value) >> 17));
case 0xe: return((value >> 13) - (value >> 16)); case 0x8e: return((value >> 16) - (value >> 13)); case 0x10e: return((value >> 13) - (value >> 16)); case 0x18e: return((value >> 16) - (value >> 13));
case 0xf: return((value >> 13) - (value >> 17)); case 0x8f: return((value >> 17) - (value >> 13)); case 0x10f: return((value >> 13) - (value >> 17)); case 0x18f: return((value >> 17) - (value >> 13));
case 0x10: return(value >> 13); case 0x90: return(-(value >> 13)); case 0x110: return(-((-value) >> 13)); case 0x190: return((-value) >> 13);
case 0x11: return((value >> 13) + (value >> 16)); case 0x91: return(-(value >> 13) - (value >> 16)); case 0x111: return(-((-value) >> 13) - ((-value) >> 16)); case 0x191: return(((-value) >> 13) + ((-value) >> 16));
case 0x12: return((value >> 13) + (value >> 15)); case 0x92: return(-(value >> 13) - (value >> 15)); case 0x112: return(-((-value) >> 13) - ((-value) >> 15)); case 0x192: return(((-value) >> 13) + ((-value) >> 15));
case 0x13: return(((value >> 12) - (value >> 14)) - (value >> 16)); case 0x93: return(((value >> 14) - (value >> 12)) + (value >> 16) ); case 0x113: return(((value >> 12) - (value >> 14)) + ((-value) >> 16) ); case 0x193: return(((value >> 14) - (value >> 12)) - ((-value) >> 16));
case 0x14: return((value >> 12) - (value >> 14)); case 0x94: return((value >> 14) - (value >> 12)); case 0x114: return((value >> 12) - (value >> 14)); case 0x194: return((value >> 14) - (value >> 12));
case 0x15: return(((value >> 12) - (value >> 14)) + (value >> 16)); case 0x95: return(((value >> 14) - (value >> 12)) + (value >> 16) ); case 0x115: return(((value >> 12) - (value >> 14)) + ((-value) >> 16) ); case 0x195: return(((value >> 14) - (value >> 12)) - ((-value) >> 16));
case 0x16: return((value >> 12) - (value >> 15)); case 0x96: return((value >> 15) - (value >> 12)); case 0x116: return((value >> 12) - (value >> 15)); case 0x196: return((value >> 15) - (value >> 12));
case 0x17: return((value >> 12) - (value >> 16)); case 0x97: return((value >> 16) - (value >> 12)); case 0x117: return((value >> 12) - (value >> 16)); case 0x197: return((value >> 16) - (value >> 12));
case 0x18: return(value >> 12); case 0x98: return(-(value >> 12)); case 0x118: return(-((-value) >> 12)); case 0x198: return((-value) >> 12);
case 0x19: return((value >> 12) + (value >> 15)); case 0x99: return(-(value >> 12) - (value >> 15)); case 0x119: return(-((-value) >> 12) - ((-value) >> 15)); case 0x199: return(((-value) >> 12) + ((-value) >> 15));
case 0x1a: return((value >> 12) + (value >> 14)); case 0x9a: return(-(value >> 12) - (value >> 14)); case 0x11a: return(-((-value) >> 12) - ((-value) >> 14)); case 0x19a: return(((-value) >> 12) + ((-value) >> 14));
case 0x1b: return(((value >> 11) - (value >> 13)) - (value >> 15)); case 0x9b: return(((value >> 13) - (value >> 11)) + (value >> 15) ); case 0x11b: return(((value >> 11) - (value >> 13)) + ((-value) >> 15) ); case 0x19b: return(((value >> 13) - (value >> 11)) - ((-value) >> 15));
case 0x1c: return((value >> 11) - (value >> 13)); case 0x9c: return((value >> 13) - (value >> 11)); case 0x11c: return((value >> 11) - (value >> 13)); case 0x19c: return((value >> 13) - (value >> 11));
case 0x1d: return(((value >> 11) - (value >> 13)) + (value >> 15)); case 0x9d: return(((value >> 13) - (value >> 11)) + (value >> 15) ); case 0x11d: return(((value >> 11) - (value >> 13)) + ((-value) >> 15) ); case 0x19d: return(((value >> 13) - (value >> 11)) - ((-value) >> 15));
case 0x1e: return((value >> 11) - (value >> 14)); case 0x9e: return((value >> 14) - (value >> 11)); case 0x11e: return((value >> 11) - (value >> 14)); case 0x19e: return((value >> 14) - (value >> 11));
case 0x1f: return((value >> 11) - (value >> 15)); case 0x9f: return((value >> 15) - (value >> 11)); case 0x11f: return((value >> 11) - (value >> 15)); case 0x19f: return((value >> 15) - (value >> 11));
case 0x20: return(value >> 11); case 0xa0: return(-(value >> 11)); case 0x120: return(-((-value) >> 11)); case 0x1a0: return((-value) >> 11);
case 0x21: return((value >> 11) + (value >> 14)); case 0xa1: return(-(value >> 11) - (value >> 14)); case 0x121: return(-((-value) >> 11) - ((-value) >> 14)); case 0x1a1: return(((-value) >> 11) + ((-value) >> 14));
case 0x22: return((value >> 11) + (value >> 13)); case 0xa2: return(-(value >> 11) - (value >> 13)); case 0x122: return(-((-value) >> 11) - ((-value) >> 13)); case 0x1a2: return(((-value) >> 11) + ((-value) >> 13));
case 0x23: return(((value >> 10) - (value >> 12)) - (value >> 14)); case 0xa3: return(((value >> 12) - (value >> 10)) + (value >> 14) ); case 0x123: return(((value >> 10) - (value >> 12)) + ((-value) >> 14) ); case 0x1a3: return(((value >> 12) - (value >> 10)) - ((-value) >> 14));
case 0x24: return((value >> 10) - (value >> 12)); case 0xa4: return((value >> 12) - (value >> 10)); case 0x124: return((value >> 10) - (value >> 12)); case 0x1a4: return((value >> 12) - (value >> 10));
case 0x25: return(((value >> 10) - (value >> 12)) + (value >> 14)); case 0xa5: return(((value >> 12) - (value >> 10)) + (value >> 14) ); case 0x125: return(((value >> 10) - (value >> 12)) + ((-value) >> 14) ); case 0x1a5: return(((value >> 12) - (value >> 10)) - ((-value) >> 14));
case 0x26: return((value >> 10) - (value >> 13)); case 0xa6: return((value >> 13) - (value >> 10)); case 0x126: return((value >> 10) - (value >> 13)); case 0x1a6: return((value >> 13) - (value >> 10));
case 0x27: return((value >> 10) - (value >> 14)); case 0xa7: return((value >> 14) - (value >> 10)); case 0x127: return((value >> 10) - (value >> 14)); case 0x1a7: return((value >> 14) - (value >> 10));
case 0x28: return(value >> 10); case 0xa8: return(-(value >> 10)); case 0x128: return(-((-value) >> 10)); case 0x1a8: return((-value) >> 10);
case 0x29: return((value >> 10) + (value >> 13)); case 0xa9: return(-(value >> 10) - (value >> 13)); case 0x129: return(-((-value) >> 10) - ((-value) >> 13)); case 0x1a9: return(((-value) >> 10) + ((-value) >> 13));
case 0x2a: return((value >> 10) + (value >> 12)); case 0xaa: return(-(value >> 10) - (value >> 12)); case 0x12a: return(-((-value) >> 10) - ((-value) >> 12)); case 0x1aa: return(((-value) >> 10) + ((-value) >> 12));
case 0x2b: return(((value >> 9) - (value >> 11)) - (value >> 13)); case 0xab: return(((value >> 11) - (value >> 9)) + (value >> 13) ); case 0x12b: return(((value >> 9) - (value >> 11)) + ((-value) >> 13) ); case 0x1ab: return(((value >> 11) - (value >> 9)) - ((-value) >> 13));
case 0x2c: return((value >> 9) - (value >> 11)); case 0xac: return((value >> 11) - (value >> 9)); case 0x12c: return((value >> 9) - (value >> 11)); case 0x1ac: return((value >> 11) - (value >> 9));
case 0x2d: return(((value >> 9) - (value >> 11)) + (value >> 13)); case 0xad: return(((value >> 11) - (value >> 9)) + (value >> 13) ); case 0x12d: return(((value >> 9) - (value >> 11)) + ((-value) >> 13) ); case 0x1ad: return(((value >> 11) - (value >> 9)) - ((-value) >> 13));
case 0x2e: return((value >> 9) - (value >> 12)); case 0xae: return((value >> 12) - (value >> 9)); case 0x12e: return((value >> 9) - (value >> 12)); case 0x1ae: return((value >> 12) - (value >> 9));
case 0x2f: return((value >> 9) - (value >> 13)); case 0xaf: return((value >> 13) - (value >> 9)); case 0x12f: return((value >> 9) - (value >> 13)); case 0x1af: return((value >> 13) - (value >> 9));
case 0x30: return(value >> 9); case 0xb0: return(-(value >> 9)); case 0x130: return(-((-value) >> 9)); case 0x1b0: return((-value) >> 9);
case 0x31: return((value >> 9) + (value >> 12)); case 0xb1: return(-(value >> 9) - (value >> 12)); case 0x131: return(-((-value) >> 9) - ((-value) >> 12)); case 0x1b1: return(((-value) >> 9) + ((-value) >> 12));
case 0x32: return((value >> 9) + (value >> 11)); case 0xb2: return(-(value >> 9) - (value >> 11)); case 0x132: return(-((-value) >> 9) - ((-value) >> 11)); case 0x1b2: return(((-value) >> 9) + ((-value) >> 11));
case 0x33: return(((value >> 8) - (value >> 10)) - (value >> 12)); case 0xb3: return(((value >> 10) - (value >> 8)) + (value >> 12) ); case 0x133: return(((value >> 8) - (value >> 10)) + ((-value) >> 12) ); case 0x1b3: return(((value >> 10) - (value >> 8)) - ((-value) >> 12));
case 0x34: return((value >> 8) - (value >> 10)); case 0xb4: return((value >> 10) - (value >> 8)); case 0x134: return((value >> 8) - (value >> 10)); case 0x1b4: return((value >> 10) - (value >> 8));
case 0x35: return(((value >> 8) - (value >> 10)) + (value >> 12)); case 0xb5: return(((value >> 10) - (value >> 8)) + (value >> 12) ); case 0x135: return(((value >> 8) - (value >> 10)) + ((-value) >> 12) ); case 0x1b5: return(((value >> 10) - (value >> 8)) - ((-value) >> 12));
case 0x36: return((value >> 8) - (value >> 11)); case 0xb6: return((value >> 11) - (value >> 8)); case 0x136: return((value >> 8) - (value >> 11)); case 0x1b6: return((value >> 11) - (value >> 8));
case 0x37: return((value >> 8) - (value >> 12)); case 0xb7: return((value >> 12) - (value >> 8)); case 0x137: return((value >> 8) - (value >> 12)); case 0x1b7: return((value >> 12) - (value >> 8));
case 0x38: return(value >> 8); case 0xb8: return(-(value >> 8)); case 0x138: return(-((-value) >> 8)); case 0x1b8: return((-value) >> 8);
case 0x39: return((value >> 8) + (value >> 11)); case 0xb9: return(-(value >> 8) - (value >> 11)); case 0x139: return(-((-value) >> 8) - ((-value) >> 11)); case 0x1b9: return(((-value) >> 8) + ((-value) >> 11));
case 0x3a: return((value >> 8) + (value >> 10)); case 0xba: return(-(value >> 8) - (value >> 10)); case 0x13a: return(-((-value) >> 8) - ((-value) >> 10)); case 0x1ba: return(((-value) >> 8) + ((-value) >> 10));
case 0x3b: return(((value >> 7) - (value >> 9)) - (value >> 11)); case 0xbb: return(((value >> 9) - (value >> 7)) + (value >> 11) ); case 0x13b: return(((value >> 7) - (value >> 9)) + ((-value) >> 11) ); case 0x1bb: return(((value >> 9) - (value >> 7)) - ((-value) >> 11));
case 0x3c: return((value >> 7) - (value >> 9)); case 0xbc: return((value >> 9) - (value >> 7)); case 0x13c: return((value >> 7) - (value >> 9)); case 0x1bc: return((value >> 9) - (value >> 7));
case 0x3d: return(((value >> 7) - (value >> 9)) + (value >> 11)); case 0xbd: return(((value >> 9) - (value >> 7)) + (value >> 11) ); case 0x13d: return(((value >> 7) - (value >> 9)) + ((-value) >> 11) ); case 0x1bd: return(((value >> 9) - (value >> 7)) - ((-value) >> 11));
case 0x3e: return((value >> 7) - (value >> 10)); case 0xbe: return((value >> 10) - (value >> 7)); case 0x13e: return((value >> 7) - (value >> 10)); case 0x1be: return((value >> 10) - (value >> 7));
case 0x3f: return((value >> 7) - (value >> 11)); case 0xbf: return((value >> 11) - (value >> 7)); case 0x13f: return((value >> 7) - (value >> 11)); case 0x1bf: return((value >> 11) - (value >> 7));
case 0x40: return(value >> 7); case 0xc0: return(-(value >> 7)); case 0x140: return(-((-value) >> 7)); case 0x1c0: return((-value) >> 7);
case 0x41: return((value >> 7) + (value >> 10)); case 0xc1: return(-(value >> 7) - (value >> 10)); case 0x141: return(-((-value) >> 7) - ((-value) >> 10)); case 0x1c1: return(((-value) >> 7) + ((-value) >> 10));
case 0x42: return((value >> 7) + (value >> 9)); case 0xc2: return(-(value >> 7) - (value >> 9)); case 0x142: return(-((-value) >> 7) - ((-value) >> 9)); case 0x1c2: return(((-value) >> 7) + ((-value) >> 9));
case 0x43: return(((value >> 6) - (value >> 8)) - (value >> 10)); case 0xc3: return(((value >> 8) - (value >> 6)) + (value >> 10) ); case 0x143: return(((value >> 6) - (value >> 8)) + ((-value) >> 10) ); case 0x1c3: return(((value >> 8) - (value >> 6)) - ((-value) >> 10));
case 0x44: return((value >> 6) - (value >> 8)); case 0xc4: return((value >> 8) - (value >> 6)); case 0x144: return((value >> 6) - (value >> 8)); case 0x1c4: return((value >> 8) - (value >> 6));
case 0x45: return(((value >> 6) - (value >> 8)) + (value >> 10)); case 0xc5: return(((value >> 8) - (value >> 6)) + (value >> 10) ); case 0x145: return(((value >> 6) - (value >> 8)) + ((-value) >> 10) ); case 0x1c5: return(((value >> 8) - (value >> 6)) - ((-value) >> 10));
case 0x46: return((value >> 6) - (value >> 9)); case 0xc6: return((value >> 9) - (value >> 6)); case 0x146: return((value >> 6) - (value >> 9)); case 0x1c6: return((value >> 9) - (value >> 6));
case 0x47: return((value >> 6) - (value >> 10)); case 0xc7: return((value >> 10) - (value >> 6)); case 0x147: return((value >> 6) - (value >> 10)); case 0x1c7: return((value >> 10) - (value >> 6));
case 0x48: return(value >> 6); case 0xc8: return(-(value >> 6)); case 0x148: return(-((-value) >> 6)); case 0x1c8: return((-value) >> 6);
case 0x49: return((value >> 6) + (value >> 9)); case 0xc9: return(-(value >> 6) - (value >> 9)); case 0x149: return(-((-value) >> 6) - ((-value) >> 9)); case 0x1c9: return(((-value) >> 6) + ((-value) >> 9));
case 0x4a: return((value >> 6) + (value >> 8)); case 0xca: return(-(value >> 6) - (value >> 8)); case 0x14a: return(-((-value) >> 6) - ((-value) >> 8)); case 0x1ca: return(((-value) >> 6) + ((-value) >> 8));
case 0x4b: return(((value >> 5) - (value >> 7)) - (value >> 9)); case 0xcb: return(((value >> 7) - (value >> 5)) + (value >> 9) ); case 0x14b: return(((value >> 5) - (value >> 7)) + ((-value) >> 9) ); case 0x1cb: return(((value >> 7) - (value >> 5)) - ((-value) >> 9));
case 0x4c: return((value >> 5) - (value >> 7)); case 0xcc: return((value >> 7) - (value >> 5)); case 0x14c: return((value >> 5) - (value >> 7)); case 0x1cc: return((value >> 7) - (value >> 5));
case 0x4d: return(((value >> 5) - (value >> 7)) + (value >> 9)); case 0xcd: return(((value >> 7) - (value >> 5)) + (value >> 9) ); case 0x14d: return(((value >> 5) - (value >> 7)) + ((-value) >> 9) ); case 0x1cd: return(((value >> 7) - (value >> 5)) - ((-value) >> 9));
case 0x4e: return((value >> 5) - (value >> 8)); case 0xce: return((value >> 8) - (value >> 5)); case 0x14e: return((value >> 5) - (value >> 8)); case 0x1ce: return((value >> 8) - (value >> 5));
case 0x4f: return((value >> 5) - (value >> 9)); case 0xcf: return((value >> 9) - (value >> 5)); case 0x14f: return((value >> 5) - (value >> 9)); case 0x1cf: return((value >> 9) - (value >> 5));
case 0x50: return(value >> 5); case 0xd0: return(-(value >> 5)); case 0x150: return(-((-value) >> 5)); case 0x1d0: return((-value) >> 5);
case 0x51: return((value >> 5) + (value >> 8)); case 0xd1: return(-(value >> 5) - (value >> 8)); case 0x151: return(-((-value) >> 5) - ((-value) >> 8)); case 0x1d1: return(((-value) >> 5) + ((-value) >> 8));
case 0x52: return((value >> 5) + (value >> 7)); case 0xd2: return(-(value >> 5) - (value >> 7)); case 0x152: return(-((-value) >> 5) - ((-value) >> 7)); case 0x1d2: return(((-value) >> 5) + ((-value) >> 7));
case 0x53: return(((value >> 4) - (value >> 6)) - (value >> 8)); case 0xd3: return(((value >> 6) - (value >> 4)) + (value >> 8) ); case 0x153: return(((value >> 4) - (value >> 6)) + ((-value) >> 8) ); case 0x1d3: return(((value >> 6) - (value >> 4)) - ((-value) >> 8));
case 0x54: return((value >> 4) - (value >> 6)); case 0xd4: return((value >> 6) - (value >> 4)); case 0x154: return((value >> 4) - (value >> 6)); case 0x1d4: return((value >> 6) - (value >> 4));
case 0x55: return(((value >> 4) - (value >> 6)) + (value >> 8)); case 0xd5: return(((value >> 6) - (value >> 4)) + (value >> 8) ); case 0x155: return(((value >> 4) - (value >> 6)) + ((-value) >> 8) ); case 0x1d5: return(((value >> 6) - (value >> 4)) - ((-value) >> 8));
case 0x56: return((value >> 4) - (value >> 7)); case 0xd6: return((value >> 7) - (value >> 4)); case 0x156: return((value >> 4) - (value >> 7)); case 0x1d6: return((value >> 7) - (value >> 4));
case 0x57: return((value >> 4) - (value >> 8)); case 0xd7: return((value >> 8) - (value >> 4)); case 0x157: return((value >> 4) - (value >> 8)); case 0x1d7: return((value >> 8) - (value >> 4));
case 0x58: return(value >> 4); case 0xd8: return(-(value >> 4)); case 0x158: return(-((-value) >> 4)); case 0x1d8: return((-value) >> 4);
case 0x59: return((value >> 4) + (value >> 7)); case 0xd9: return(-(value >> 4) - (value >> 7)); case 0x159: return(-((-value) >> 4) - ((-value) >> 7)); case 0x1d9: return(((-value) >> 4) + ((-value) >> 7));
case 0x5a: return((value >> 4) + (value >> 6)); case 0xda: return(-(value >> 4) - (value >> 6)); case 0x15a: return(-((-value) >> 4) - ((-value) >> 6)); case 0x1da: return(((-value) >> 4) + ((-value) >> 6));
case 0x5b: return(((value >> 3) - (value >> 5)) - (value >> 7)); case 0xdb: return(((value >> 5) - (value >> 3)) + (value >> 7) ); case 0x15b: return(((value >> 3) - (value >> 5)) + ((-value) >> 7) ); case 0x1db: return(((value >> 5) - (value >> 3)) - ((-value) >> 7));
case 0x5c: return((value >> 3) - (value >> 5)); case 0xdc: return((value >> 5) - (value >> 3)); case 0x15c: return((value >> 3) - (value >> 5)); case 0x1dc: return((value >> 5) - (value >> 3));
case 0x5d: return(((value >> 3) - (value >> 5)) + (value >> 7)); case 0xdd: return(((value >> 5) - (value >> 3)) + (value >> 7) ); case 0x15d: return(((value >> 3) - (value >> 5)) + ((-value) >> 7) ); case 0x1dd: return(((value >> 5) - (value >> 3)) - ((-value) >> 7));
case 0x5e: return((value >> 3) - (value >> 6)); case 0xde: return((value >> 6) - (value >> 3)); case 0x15e: return((value >> 3) - (value >> 6)); case 0x1de: return((value >> 6) - (value >> 3));
case 0x5f: return((value >> 3) - (value >> 7)); case 0xdf: return((value >> 7) - (value >> 3)); case 0x15f: return((value >> 3) - (value >> 7)); case 0x1df: return((value >> 7) - (value >> 3));
case 0x60: return(value >> 3); case 0xe0: return(-(value >> 3)); case 0x160: return(-((-value) >> 3)); case 0x1e0: return((-value) >> 3);
case 0x61: return((value >> 3) + (value >> 6)); case 0xe1: return(-(value >> 3) - (value >> 6)); case 0x161: return(-((-value) >> 3) - ((-value) >> 6)); case 0x1e1: return(((-value) >> 3) + ((-value) >> 6));
case 0x62: return((value >> 3) + (value >> 5)); case 0xe2: return(-(value >> 3) - (value >> 5)); case 0x162: return(-((-value) >> 3) - ((-value) >> 5)); case 0x1e2: return(((-value) >> 3) + ((-value) >> 5));
case 0x63: return(((value >> 2) - (value >> 4)) - (value >> 6)); case 0xe3: return(((value >> 4) - (value >> 2)) + (value >> 6) ); case 0x163: return(((value >> 2) - (value >> 4)) + ((-value) >> 6) ); case 0x1e3: return(((value >> 4) - (value >> 2)) - ((-value) >> 6));
case 0x64: return((value >> 2) - (value >> 4)); case 0xe4: return((value >> 4) - (value >> 2)); case 0x164: return((value >> 2) - (value >> 4)); case 0x1e4: return((value >> 4) - (value >> 2));
case 0x65: return(((value >> 2) - (value >> 4)) + (value >> 6)); case 0xe5: return(((value >> 4) - (value >> 2)) + (value >> 6) ); case 0x165: return(((value >> 2) - (value >> 4)) + ((-value) >> 6) ); case 0x1e5: return(((value >> 4) - (value >> 2)) - ((-value) >> 6));
case 0x66: return((value >> 2) - (value >> 5)); case 0xe6: return((value >> 5) - (value >> 2)); case 0x166: return((value >> 2) - (value >> 5)); case 0x1e6: return((value >> 5) - (value >> 2));
case 0x67: return((value >> 2) - (value >> 6)); case 0xe7: return((value >> 6) - (value >> 2)); case 0x167: return((value >> 2) - (value >> 6)); case 0x1e7: return((value >> 6) - (value >> 2));
case 0x68: return(value >> 2); case 0xe8: return(-(value >> 2)); case 0x168: return(-((-value) >> 2)); case 0x1e8: return((-value) >> 2);
case 0x69: return((value >> 2) + (value >> 5)); case 0xe9: return(-(value >> 2) - (value >> 5)); case 0x169: return(-((-value) >> 2) - ((-value) >> 5)); case 0x1e9: return(((-value) >> 2) + ((-value) >> 5));
case 0x6a: return((value >> 2) + (value >> 4)); case 0xea: return(-(value >> 2) - (value >> 4)); case 0x16a: return(-((-value) >> 2) - ((-value) >> 4)); case 0x1ea: return(((-value) >> 2) + ((-value) >> 4));
case 0x6b: return(((value >> 1) - (value >> 3)) - (value >> 5)); case 0xeb: return(((value >> 3) - (value >> 1)) + (value >> 5) ); case 0x16b: return(((value >> 1) - (value >> 3)) + ((-value) >> 5) ); case 0x1eb: return(((value >> 3) - (value >> 1)) - ((-value) >> 5));
case 0x6c: return((value >> 1) - (value >> 3)); case 0xec: return((value >> 3) - (value >> 1)); case 0x16c: return((value >> 1) - (value >> 3)); case 0x1ec: return((value >> 3) - (value >> 1));
case 0x6d: return(((value >> 1) - (value >> 3)) + (value >> 5)); case 0xed: return(((value >> 3) - (value >> 1)) + (value >> 5) ); case 0x16d: return(((value >> 1) - (value >> 3)) + ((-value) >> 5) ); case 0x1ed: return(((value >> 3) - (value >> 1)) - ((-value) >> 5));
case 0x6e: return((value >> 1) - (value >> 4)); case 0xee: return((value >> 4) - (value >> 1)); case 0x16e: return((value >> 1) - (value >> 4)); case 0x1ee: return((value >> 4) - (value >> 1));
case 0x6f: return((value >> 1) - (value >> 5)); case 0xef: return((value >> 5) - (value >> 1)); case 0x16f: return((value >> 1) - (value >> 5)); case 0x1ef: return((value >> 5) - (value >> 1));
case 0x70: return(value >> 1); case 0xf0: return(-(value >> 1)); case 0x170: return(-((-value) >> 1)); case 0x1f0: return((-value) >> 1);
case 0x71: return((value >> 1) + (value >> 4)); case 0xf1: return(-(value >> 1) - (value >> 4)); case 0x171: return(-((-value) >> 1) - ((-value) >> 4)); case 0x1f1: return(((-value) >> 1) + ((-value) >> 4));
case 0x72: return((value >> 1) + (value >> 3)); case 0xf2: return(-(value >> 1) - (value >> 3)); case 0x172: return(-((-value) >> 1) - ((-value) >> 3)); case 0x1f2: return(((-value) >> 1) + ((-value) >> 3));
case 0x73: return(((value >> 0) - (value >> 2)) - (value >> 4)); case 0xf3: return(((value >> 2) - (value >> 0)) + (value >> 4) ); case 0x173: return(((value >> 0) - (value >> 2)) + ((-value) >> 4) ); case 0x1f3: return(((value >> 2) - (value >> 0)) - ((-value) >> 4));
case 0x74: return((value >> 0) - (value >> 2)); case 0xf4: return((value >> 2) - (value >> 0)); case 0x174: return((value >> 0) - (value >> 2)); case 0x1f4: return((value >> 2) - (value >> 0));
case 0x75: return(((value >> 0) - (value >> 2)) + (value >> 4)); case 0xf5: return(((value >> 2) - (value >> 0)) + (value >> 4) ); case 0x175: return(((value >> 0) - (value >> 2)) + ((-value) >> 4) ); case 0x1f5: return(((value >> 2) - (value >> 0)) - ((-value) >> 4));
case 0x76: return((value >> 0) - (value >> 3)); case 0xf6: return((value >> 3) - (value >> 0)); case 0x176: return((value >> 0) - (value >> 3)); case 0x1f6: return((value >> 3) - (value >> 0));
case 0x77: return((value >> 0) - (value >> 4)); case 0xf7: return((value >> 4) - (value >> 0)); case 0x177: return((value >> 0) - (value >> 4)); case 0x1f7: return((value >> 4) - (value >> 0));
case 0x78: return(value >> 0); case 0xf8: return(-(value >> 0)); case 0x178: return(-((-value) >> 0)); case 0x1f8: return((-value) >> 0);
case 0x79: return((value >> 0) + (value >> 3)); case 0xf9: return(-(value >> 0) - (value >> 3)); case 0x179: return(-((-value) >> 0) - ((-value) >> 3)); case 0x1f9: return(((-value) >> 0) + ((-value) >> 3));
case 0x7a: return((value >> 0) + (value >> 2)); case 0xfa: return(-(value >> 0) - (value >> 2)); case 0x17a: return(-((-value) >> 0) - ((-value) >> 2)); case 0x1fa: return(((-value) >> 0) + ((-value) >> 2));
case 0x7b: return(((value << 1) - (value >> 1)) - (value >> 3)); case 0xfb: return(((value >> 1) - (value << 1)) + (value >> 3) ); case 0x17b: return(((value << 1) - (value >> 1)) + ((-value) >> 3) ); case 0x1fb: return(((value >> 1) - (value << 1)) - ((-value) >> 3));
case 0x7c: return((value << 1) - (value >> 1)); case 0xfc: return((value >> 1) - (value << 1)); case 0x17c: return((value << 1) - (value >> 1)); case 0x1fc: return((value >> 1) - (value << 1));
case 0x7d: return(((value << 1) - (value >> 1)) + (value >> 3)); case 0xfd: return(((value >> 1) - (value << 1)) + (value >> 3) ); case 0x17d: return(((value << 1) - (value >> 1)) + ((-value) >> 3) ); case 0x1fd: return(((value >> 1) - (value << 1)) - ((-value) >> 3));
case 0x7e: return((value << 1) - (value >> 2)); case 0xfe: return((value >> 2) - (value << 1)); case 0x17e: return((value << 1) - (value >> 2)); case 0x1fe: return((value >> 2) - (value << 1));
case 0x7f: return((value << 1) - (value >> 3)); case 0xff: return((value >> 3) - (value << 1)); case 0x17f: return((value << 1) - (value >> 3)); case 0x1ff: return((value >> 3) - (value << 1));
		default: return(0);
		}

	}*/

abstractMesh::abstractMesh(): nbpts(0){}

void abstractMesh::insertTriangle(unsigned short a, unsigned short b, unsigned short c){
	abstractMeshFace newface;
	newface.points[0] =a ;newface.points[1] =b ;newface.points[2] =c;
	memset(newface.adjacent,'\xFF',sizeof(unsigned short)*3);
	memset(&(newface.tag),'\xFF',sizeof(int));
	mesh.push_back(newface);
}

void abstractMesh::insertNew2dGrid(unsigned short x, unsigned short y){
	int i,j;
	abstractMeshFace newface;
	memset(newface.adjacent,'\xFF',sizeof(unsigned short)*3);
	memset(&(newface.tag),'\xFF',sizeof(int));
	for(j=0;j<y;j++){
		for(i=0;i<x;i++){
			newface.points[0] = i + j* (x+1);
			newface.points[1] = (i+1)+ j* (x+1);
			newface.points[2] = i+ (j+1)* (x+1);
			mesh.push_back(newface);
			newface.points[0] = newface.points[2];
			newface.points[2] = (i+1)+ (j+1)* (x+1);
			mesh.push_back(newface);
		}
	}
	nbpts += (x+1) * (y+1);
}
/*
void abstractMesh::insertNew2dGridHole(unsigned short x, unsigned short y, unsigned short hole){
	bool safe;
	for(j=0;j<y;j++){
		safe = (j - (y>>1)) >
		for(i=0;i<x;i++){
			newface.points[0] = i + j* (x+1);
			newface.points[1] = (i+1)+ j* (x+1);
			newface.points[2] = i+ (j+1)* (x+1);
			mesh.push_back(newface);
			newface.points[0] = newface.points[2];
			newface.points[2] = (i+1)+ (j+1)* (x+1);
			mesh.push_back(newface);
		}
	}

	nbpts += (x+1) * (y+1) - (hole-1) * (hole-1);

}*/

void abstractMesh::fillajacents(){

}

int OwnerRegisterCmp::operator()(void*  const a , void* const b) const{
	return( (a > b) ? 1 : (a == b ? 0 : -1));
	}

OwnerRegisterEntry::OwnerRegisterEntry(void *owner, void* start, void* end):  key(start), endkey(end), target(owner){

	}
void* OwnerRegisterEntry::getIndex(){return(key);}



/*
static void loadBMP( char const * const path, DataGrid<Tuple<unsigned char, 3>,2> &im){

	//printf("Opening %s\n",path); fflush(stdout);
	FILE* f = fopen(path,"rb+"); if (f == NULL) exit(123);
	char buffer[65536];
	fread(buffer,sizeof(char),54,f);
	uint32_t dims[2];
	dims[0] = *((int*)(buffer + 18));
	dims[1] = *((int*)(buffer + 22));

	unsigned short format = *((unsigned short*)(buffer + 28));

	im.setSizes(dims);
	uint32_t i,j;
	if (format == 24){

	//printf("size: %i,%i\n", im.dims[0], im.dims[1]);fflush(stdout);
	for(i=0;i< dims[1];i++){
		fread( &(im[0])+ i*dims[0] ,sizeof(char),dims[0]*3,f);
		if ((dims[0]*3) % 4 != 0) fread(buffer,sizeof(char),4-((dims[0]*3)%4),f);
	}
	}else{
		for(i=0;i< dims[1];i++) for(j=0;j< dims[0];j++){
			fread( &(im[0]) + i*dims[0] + j ,sizeof(char),3,f);
			fread(buffer,sizeof(char), (format / 8) -3,f);
			if ((dims[0]*(format / 8)) % 4 != 0) fread(buffer,sizeof(char),4-((dims[0]*(format / 8))%4),f);
		}

	}


	// printf("Opened %s sucessfully\n",path);

	fclose(f);



	}
static void saveBMP(const char* path, const DataGrid<Tuple<unsigned char, 3>,2> &im){
		FILE* f = fopen(path,"wb+");

				char header[54];
	char shead[] = {'B', 'M',
					0,0,0,0,
					'L', 'F', 'H', 'F',
					54,0,0,0,
					40,0,0,0,
					0,0,0,0,
					0,0,0,0,
					1,0,
					24,0, // image depth
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0};
	memcpy(header,shead, sizeof(char)*54);
	int* write;

	write = (int*)(header + 18); (*write) = im.dims[0];
	write = (int*)(header + 22); (*write) = im.dims[1];
	write = (int*)(header + 2); (*write) = 54 + im.dims[0] * im.dims[1] * 3;
	uint32_t i;
	fwrite(header,sizeof(char),54,f);
	memset(header,'\0',sizeof(char)*4);

	for(i=0;i< im.dims[1];i++){
		fwrite(im.data + i*im.dims[0] ,sizeof(char),im.dims[0]*3,f);
		if ((im.dims[0]*3) % 4 != 0) fwrite(header,sizeof(char),4-((im.dims[0]*3)%4),f);
	}
	fclose(f);


	}*/

void FileIO::loadBMP3(const char* path, DataGrid<unsigned char,3> &im){

	//printf("Opening %s\n",path); fflush(stdout);
	FILE* f = fopen(path,"rb+"); if (f == NULL) exit(123);
	unsigned char buffer[65536];
	if (fread(buffer,sizeof(char),54,f) != 54) return;
	uint32_t dims[3];
	dims[1] = *((int*)(buffer + 18));
	dims[2] = *((int*)(buffer + 22));

	unsigned short format = *((unsigned short*)(buffer + 28));
	dims[0] = format / 8;

	im.setSizes(dims);
	uint32_t i,j;
	int32_t ret;
	if (format == 24){
	//printf("size: %i,%i\n", im.dims[0], im.dims[1]);fflush(stdout);
	for(i=0;i< dims[2];i++){
		ret = fread( &(im.data[0])+ i*dims[1]*dims[0] ,sizeof(char),dims[1]*dims[0],f);
		if ((dims[1]*dims[0]) % 4 != 0) ret = fread(buffer,sizeof(char),4-((dims[1]*dims[0])%4),f);
	}
	}else{
		for(i=0;i< dims[2];i++) for(j=0;j< dims[1];j++){
			ret = fread( &(im.data[0]) + i*dims[1]*dims[0] + j ,sizeof(char),dims[0],f);
			if (dims[0] != 3) ret = fread(buffer,sizeof(char), dims[0] -3,f);
			if ((dims[1]*dims[0]) % 4 != 0) ret = fread(buffer,sizeof(char),4-((dims[1]*dims[0])%4),f);
		}

	}
	// printf("Opened %s sucessfully\n",path);
	printf("ret is %i\n", ret); // no more warning from this obsolete code!
	fclose(f);
}
void FileIO::saveBMP3(const char* const path,const DataGrid<unsigned char, 3> &im){
    FILE* f = fopen(path,"wb+");
    char header[54];
	char shead[] = {'B', 'M',
					0,0,0,0,
					'L', 'F', 'H', 'F',
					54,0,0,0,
					40,0,0,0,
					0,0,0,0,
					0,0,0,0,
					1,0,
					24,0, // image depth
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0};
	memcpy(header,shead, sizeof(char)*54);
	int* write;

	write = (int*)(header + 18); (*write) = im.dims[1];
	write = (int*)(header + 22); (*write) = im.dims[2];
	write = (int*)(header + 2); (*write) = 54 + im.dims[1] * im.dims[2] * 3;
	uint32_t i;
	fwrite(header,sizeof(char),54,f);
	memset(header,'\0',sizeof(char)*4);

	for(i=0;i< im.dims[2];i++){
		fwrite(im.data + i*im.dims[1] ,sizeof(char),im.dims[1]*3,f);
		if ((im.dims[1]*3) % 4 != 0) fwrite(header,sizeof(char),4-((im.dims[1]*3)%4),f);
	}
	fclose(f);


	}



	/*
void EventMaintainer::ThreadProgram(){
   while(true){
    sync_event_lock++;
    if (sync_event_lock == 1){
        // unique thread
        // manages the unique events and the locks
        while(true){
      //      printf("is main thread!\n");
//            SwitchToThread();
        }
    }
    sync_event_lock--;
    if (sync_event_lock != 0) break;
    }
        while(true){
     //       printf("is local thread!\n");
     //      SwitchToThread();
        }
}
*/

template< >	double QuadraticBound_old<double>::minRoot() const{
	//P(x) = (S+E)/2 + x(E-S)/2 +- K (1 - x2)
/*	if (en * st < 0){ // root in range guarantied

	}else{ // two roots in range or none

	}

	double c = st + en; // /4
	double a = err*2.0f;
	double d = b*b + (c+a) * a;
	if (d >= 0.0f){

	}

	d = b*b - (c-a) * a;
	if (d >= 0.0f){

	}*/
	if (err == 0.0f){
		if (st * en > 0) return(1);
		else return((st + en) / (st - en));
	}
	double c = (st + en);
	double b = st - en; // /-2
	if (st < 0){ // care about a > 0 too
		c *= -1.0f;
		b *= -0.5f;
	}else b *= 0.5f;
	//printf("%f\t%f\n",c,b);

	double e = err *2.0f;

	c = b*b - e * (c-e);
	//printf("%f\t%f\n",e,c);
	if (c < 0) return(1);
	else return( (b - sqrt(c) ) / e);
}

template< >	double QuarticBound<double>::minRoot() const{
	double grec[3];
	double tmp;
	grec[2] = (data[0] > 0) ? 0.5f / (data[4] - err) : 0.5f / (data[2] + err);
	tmp = grec[2] * data[3];

	grec[0] = ((((-3.0f / 16.0f)* tmp * data[3] + 0.5f*data[2])*tmp - data[1]) * tmp + 2.0f * data[0]) * grec[2];
	grec[1] = ((tmp*data[3] -0.5f*data[2])*tmp + 2.0f*data[1])* grec[2];

	grec[2] = ((-3.0f/2.0f)*data[3]*data[3]*grec[2] + 2.0f * data[2])/grec[2];

	tmp = (grec[2]*grec[2] / -12.0f - grec[0])/ 3.0f ; // P / 3
	grec[0] = (((grec[2] - (grec[0]*grec[0]/16.0f)) * grec[0] / 3.0f) - grec[1]*grec[1]/8.0f)/2.0f; // Q / 2
	double discr = grec[0] * grec[0] + tmp*tmp*tmp;

    return(discr);
}



/*

    AliasHost::~AliasHost(){
       if (object_alias){
            uint32_t ite = static_alias_bank.pointers.find(object_alias);
            if (static_alias_bank.pointers.deref(ite).first == (void*) this) {
                if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");LFH_exit(1);}
                static_alias_bank.pointers.deref(ite).first = NULL;
                if (static_alias_bank.pointers.deref(ite).second == 0) {fprintf(stderr,"Alias entry was not deleted as expected!\n");LFH_exit(1);}
            } // else entry was stolen by another Host
        }
    }


    AliasHost::AliasHost(const AliasHost& other): object_alias(other.object_alias){ // steals entry!
        if (object_alias){
            uint32_t ite = static_alias_bank.pointers.find(object_alias);
            if (static_alias_bank.pointers.deref(ite).first == (void*) &other) static_alias_bank.pointers.deref(ite).first = (void*) this;
        }
    }

    AliasHost& AliasHost::operator=(const AliasHost& other){
        if (object_alias){
            uint32_t ite = static_alias_bank.pointers.find(object_alias);
            if (static_alias_bank.pointers.deref(ite).first == (void*) this) {
                if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");LFH_exit(1);}
                static_alias_bank.pointers.deref(ite).first = NULL;
                if (static_alias_bank.pointers.deref(ite).second == 0) {fprintf(stderr,"Alias entry was not deleted as expected!\n");LFH_exit(1);}
            } // else entry was stolen by another Host
        }
        object_alias = other.object_alias;
        if (object_alias){
            uint32_t ite = static_alias_bank.pointers.find(object_alias);
            if (static_alias_bank.pointers.deref(ite).first == (void*) &other) static_alias_bank.pointers.deref(ite).first = (void*) this;
        }
        return(*this);
    }
*/
/*
Alias::~Alias(){
    if (alias){
        uint32_t ite = central_alias_bank.find(alias);
        if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");LFH_exit(1);}
        central_alias_bank.deref(ite).first = NULL;
    }
}*/

AliasHost::AliasHost(bool autoAssign){if (autoAssign) this->setAHAlias();}

void AliasHost::setAHAlias(uint32_t r){
    uint32_t ite;
    if (r == 0){
        do{
        ExOp::toRand(r);
        while(r == 0) ExOp::toRand(r);
        ite = AliasBank.find(r);
        }while(ite != 0xFFFFFFFF);
    }
    AliasBank[r] = pair<void*, uint32_t>(this,0u);
    AliasOf[this] = r;
}

AliasHost::~AliasHost(){
    uint32_t alias = AliasOf[this];
    if (alias){
    AliasBank[alias].first = NULL;
    AliasOf.erase(this);
    }
}



void FIFOQueue<void>::runAll(){
    while(semaphore != async_read){
        (*buffer[async_read])();
        delete(buffer[async_read]);
        async_read++;
        if (async_read & 0x100){
            async_read &= 0xFF;
            semaphore -= 0x100;
        }
    }
}
			template<class C, int nbc>
	void imageIO<C,nbc>::importBMP(char* path, DataGrid<Tuple<C, nbc>,2> &out){

		}

	template<class C, int nbc>
	void imageIO<C,nbc>::exportBMP(char* path, DataGrid<Tuple<C, nbc>,2> &in){
			char header[54];
	char shead[] = {'B', 'M',
					0,0,0,0,
					'L', 'F', 'H', 'F',
					54,0,0,0,
					40,0,0,0,
					0,0,0,0,
					0,0,0,0,
					1,0,
					24,0, // image depth
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0,
					0,0,0,0};
	 // "BM    LFHF\0\0\0\54\0\0\0(        \255\255\                                "
	memcpy(header,shead, sizeof(char)*54);
	int* write;
	write = (int*)(header + 18); (*write) = in.dims[0];
	write = (int*)(header + 22); (*write) = in.dims[1];
	write = (int*)(header + 2); (*write) = 54 + in.dims[0] * in.dims[1] * 3;
	FILE* g = fopen(path,"wb+");

	fwrite(header,sizeof(char),54,g);

	int i;


	memset(header,'\0',sizeof(char)*4);
//	for(i=0;i< in.dims[1];i++){
//		fwrite(&data[i*sizex*3],sizeof(char),sizex*channels,g);
//		if ((sizex*3) % 4 != 0) fwrite(header,sizeof(char),4-((sizex*3)%4),g);
//	}
	fclose(g);

		}


void RunningStatistics<true>::begin(const char* label){
    string s(label);
    uint32_t ite = procmap.find(s);
    if (ite == 0xFFFFFFFF){
        procmap[s] = runs.getSize();
        runs.push_back();
        ExOp::toZero(runs[runs.getSize()-1u]);
        ite = procmap.find(s);
        procnames.push_back(s);
    }
    scope.push_back(KeyElem<uint32_t,uint32_t>(clock(), procmap.deref(ite)));
    }
void RunningStatistics<true>::end(const char* label){
    string s(label);
    uint32_t ite = procmap.find(s);
    if (ite == 0xFFFFFFFF){fprintf(stderr, "Debug error: process %s ended before it started! (wrong name?)\n", label); LFH_exit(1);}
    ite = procmap.deref(ite);
    if (scope[scope.getSize()-1u].d != ite){
        {fprintf(stderr, "Debug error: process %s ended within the scope of %s\n", label, procnames[scope[scope.getSize()-1u].d].c_str()); LFH_exit(1);}
    }
    runs[ite] += WeightElem<double,2>(clock() - scope[scope.getSize()-1u].k);
    scope.pop_back();
    }
void RunningStatistics<true>::report(bool exit_with_error){
        printf("\n\nRunning Statistics Report:\n");
        uint32_t i;
        for(i=0;i<runs.getSize();i++){
            printf("Process %s:\n", procnames[i].c_str());
            printf("\tRan %i times\n", (uint32_t)runs[i].w[0]);
            printf("\taverage milli %f +- %f\n", runs[i].getMean(), sqrt(runs[i].getVar()));
        }
        if (exit_with_error) LFH_exit(1);
    }

ScriptHead::ScriptHead(){
    parse_dictionary[string("if")] = &type; //NULL;
    parse_dictionary[string("else")] = &type; //NULL;
    parse_dictionary[string("unsigned")] = &type; //NULL;
    parse_dictionary[string("char")] = &type;
    parse_dictionary[string("int")] = &type;
    parse_dictionary[string("short")] = &type;
    parse_dictionary[string("double")] = &type;
    parse_dictionary[string("float")] = &type;
    }

Vector<Instructions> ScriptHead::parseScript(const char* script){ Vector<Instructions> fout;
    const char* cscr = script;
    stack<const char*> command;
    stack<uint32_t> operat;
    const char* escr;
    string cs;
    uint32_t j;
    while(*cscr != '\0'){
        while((*cscr == ' ')||(*cscr == '\n')||(*cscr == '\t')) cscr++;
        if (*cscr == '\0') break;
        if ((*cscr >= 'A')&&(*cscr <= 'z')&&((*cscr <='Z')|| (*cscr >='a')|| (*cscr == '_'))){ // function or variable
            escr = cscr++;
            while(((*cscr >= 'A')&&(*cscr <= 'z')&&((*cscr <='Z')|| (*cscr >='a')|| (*cscr =='_')))||((*cscr >= '0')&&(*cscr <= '9'))) cscr++;
                cs = string(escr, (uint32_t) (cscr - escr));
                j = parse_dictionary.find(cs);
                if (j != 0xFFFFFFFF){
                    printf(" found %s in dictionary\n", cs.c_str());
                }
            } else if ((*cscr >= '*')&&(*cscr <= '/')) {
                if (*cscr == '+') operat.push(1);
                if (*cscr == '-') operat.push(2);
                if (*cscr == '*') operat.push(3);
                if (*cscr == '/') operat.push(4);
                // , and . here too!
                cscr++;
            } else cscr++;
        }
    return fout;
    }

    // function calls
    // assignents
    // branching
    void ScriptHead::runScript( Vector<Instructions> &code){
        uint32_t cur = 0;
        void* buffer[16];
        uint32_t i,j;
        while(true){
            switch(code[cur++].v){
                case 0: // function call
                    i = (code[cur++].v);
                    for(j=0;j<i;j++) buffer[j] = code[++cur].p;
                    (*(ScriptNode *)code[++cur].p)(buffer);
                break;
            }
        }
    }

    void GaussianRegression::EMinit(uint32_t size, uint32_t nbstates){
        means.setSize(size);
        stds.setSize(size);
        for(uint32_t i=0;i<size;i++) {stds[i].setSize(nbstates); means[i].setSize(nbstates);}
    }
    void GaussianRegression::EMregist(const Tuple<double, 0u>&){
        // expect step (haha)
        Tuple<double,0u> wei; wei.setSize(means[0].getSize());
        ExOp::toRand(wei); wei *= wei;

        // stats for M-step
       // uint32_t s;
  //      for(s=0;s<means.size();s++){


   //     }

    }
void GaussianRegression::EMfinit(){

}

void GaussianProcessMK2::checkDomain(CurvatureSearchScope &learnscope){
	if (data[2] <= 0.0) {
	//	printf("signal got negative! %e", data[2]);
		learnscope.makeGuessStrictlyPositive(data[2], 2);
//		printf("-> %e\n", data[2]);
	}
	if (data[1] <= 0.0){
	//	printf("noise got negative!\n");
		learnscope.makeGuessStrictlyPositive(data[1],1);
	}
	unsigned int i,k;
	for(k=0,i=0;i<scale.getSize();i++,k++){
		k += i;
		if (scale.data[k] <= 0.0) learnscope.makeGuessStrictlyPositive(scale.data[k], 3+k);
	}
}

void GaussianProcessMK2::checkDomain(AdaptiveLearningScope &learnscope){
	if (learnscope.makeGuessStrictlyPositive(data[2], 2)) {
		printf("signal got negative! fixed to %e", data[2]);
	}
	if (learnscope.makeGuessStrictlyPositive(data[1], 1)){
		printf("noise got negative! fixed to %e", data[1]);
	}
	unsigned int i,k;
	for(k=0,i=0;i<scale.getSize();i++,k++){
		k += i;
		if (learnscope.makeGuessStrictlyPositive(scale.data[k], 3+k)) {
			printf("diago got negative! fixed to %e\n", scale.data[k]);
		}
	}
}

double GaussianProcessMK2::addLLTimeDerivative(Tuple<double>& fout, RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, bool verbose) const{
	unsigned int k,m;
	Tuple<double> dif = this->setKernel3(observ, timestamp);

	short posdefclass;
	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	LL += cached_kernel.log_determinant(&posdefclass);
	if (posdefclass != 1) return sqrt(-1.0);
	double altLL = 0.0;

	if (verbose){
		Tuple<double> andback = cached_kernel.leftDivision(dif);
		for(unsigned int i =0;i<dif.getSize();i++) altLL += andback[i] * dif[i];
		printf("LL %e vs %e  %e error\n", LL, altLL, (LL - altLL) /LL);
		andback = cached_kernel * andback;
		altLL = 0.0;
		for(unsigned int i =0;i<dif.getSize();i++) altLL += fabs((andback[i] - dif[i]) / dif[i]);
		printf("divi error %e\n", altLL);
		altLL = 0.0;
	}

	LL *= -0.5;
	dif = cached_kernel.leftDivision(dif);

	Tuple< Tuple<double> > ders; ders.setSize(cached_kernel.getSize() * scale.getSize());
	Tuple<double> tmp; tmp.setSize(scale.getSize());

	for(m=0;m<ders.getSize();m++) ders[m].setSize(cached_kernel.getSize()).toZero();
	Tuple<uint32_t, 2u> coor, coor2;
	for(k=0,coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
		for(coor2[1]=0;coor2[1]<coor[1];coor2[1]++,k++){ // i and j are sampleIDs
			for(coor[0]=0;coor[0]<scale.getSize();coor[0]++) {coor2[0] = coor[0]; tmp[coor[0]] = timestamp[coor2] - timestamp[coor];}
			//printf("pre\n");tmp.show();
			tmp = scale * tmp; // projection
			//printf("post\n");tmp.show();
			for(m=0;m<scale.getSize();m++){ // n and m are dimensions
				ders[coor[1] * scale.getSize() + m][coor2[1]] -= cached_kernel.data[k] * tmp[m];
				ders[coor2[1] * scale.getSize() + m][coor[1]] += cached_kernel.data[k] * tmp[m];
			}
		}
		k++;
	}

	for(k=0;k<ders.getSize();k++) {
		altLL =  dif[0] * ders[k][0];
		for(m=1;m<cached_kernel.getSize();m++) altLL += dif[m] * ders[k][m];
		altLL *= dif[(k / scale.getSize())];

		ders[k] = cached_kernel.leftDivision(ders[k]);
		altLL = ders[k][(k / scale.getSize())] - altLL;
		fout[k] += 2.0 * altLL;
	}
return LL;}
/*
double GaussianProcessMK2::getMeanAndDeriv_p1(double extra_state, double weight, Tuple<double>& statederiv, RemoteMemory<const double,2u> delta_state, double* buffer)const{
	// writes nbdim + 2 derivatives
	if (delta_state.dims[1] != cached_projection.getSize()) {printf("Cached projection is missing!!! %i, %i\n", delta_state.dims[1], cached_projection.getSize()); return 0.0;}
	double fout = 0.0;
	double current;
	int i =0;
	current = exp(-scale.Xformed_inner_product(delta_state.selectSlice(i,1)));
	current *= cached_projection[i];
	fout = current;
	Tuple<double> res = scale * delta_state.selectSlice(i,1);
	if (res.getSize() !=scale.getSize()) {printf("multiplication fail!! %i, %i\n", res.getSize(), scale.getSize()); return 0.0;}
	for(int j=0;j< scale.getSize();j++) buffer[j] = res[j] * current;
	for(i++;i< delta_state.dims[1];i++) {
		current = exp(-scale.Xformed_inner_product(delta_state.selectSlice(i,1)));
		current *= cached_projection[i];
		fout += current;
		res = scale * delta_state.selectSlice(i,1);
		for(int j=0;j< scale.getSize();j++) buffer[j] += res[j] * current;
	}
	current = exp(-extra_state*extra_state);
	fout = weight * exp(this->getMean() + this->getSignal() * fout * current);
	current *= fout * this->getSignal() * 2.0;
	for(i=0;i< scale.getSize();i++) statederiv[i] -= current * buffer[i]; // exp
	statederiv[i] -= current * extra_state;
return fout;}*/

double GaussianProcessMK2::getMeanAndDeriv_p1(double extra_state, double weight, Tuple<double>& statederiv, RemoteMemory<const double,2u> delta_state, double* buffer)const{
	// writes nbdim + 2 derivatives
	if (delta_state.dims[1] != cached_projection.getSize()) {printf("Cached projection is missing!!! %i, %i\n", delta_state.dims[1], cached_projection.getSize()); return 0.0;}
	double fout = 0.0;
	double current;
	uint32_t i =0;
	current = exp(-extra_state * extra_state - scale.Xformed_inner_product(delta_state.selectSlice(i,1)));
	current *= cached_projection[i];
	fout = current;
	double exbuf = -2.0 * current * extra_state;
	Tuple<double> res = scale * delta_state.selectSlice(i,1);
	if (res.getSize() !=scale.getSize()) {printf("multiplication fail!! %i, %i\n", res.getSize(), scale.getSize()); return 0.0;}
	for(uint32_t j=0;j< scale.getSize();j++) buffer[j] = res[j] * current;
	for(i++;i< delta_state.dims[1];i++) {
		current = exp(-extra_state * extra_state - scale.Xformed_inner_product(delta_state.selectSlice(i,1)));
		current *= cached_projection[i];
		fout += current;
		res = scale * delta_state.selectSlice(i,1);
		exbuf -= 2.0 * current * extra_state;
		for(uint32_t j=0;j< scale.getSize();j++) buffer[j] += res[j] * current;
	}
	fout = weight * exp(this->getMean() + this->getSignal() * fout);
	current = fout * this->getSignal();
	for(i=0;i< scale.getSize();i++) statederiv[i] -= 2.0 * current * buffer[i]; // exp
	statederiv[i] += current * exbuf;
return fout;}
/*
double GaussianProcessMK2::addLLTimeDerivative2(Tuple<double>& fout, double value, Tuple<double> state, RemoteMemory<const double,2u> timestamp, bool verbose) const{
	unsigned int k,m,n;

	short posdefclass;
	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	LL += cached_kernel.log_determinant(&posdefclass);
	if (posdefclass != 1) return sqrt(-1.0);
	double altLL = 0.0;

	if (verbose){
		Tuple<double> andback = cached_kernel.leftDivision(dif);
		for(unsigned int i =0;i<dif.getSize();i++) altLL += andback[i] * dif[i];
		printf("LL %e vs %e  %e error\n", LL, altLL, (LL - altLL) /LL);
		andback = cached_kernel * andback;
		altLL = 0.0;
		for(unsigned int i =0;i<dif.getSize();i++) altLL += fabs((andback[i] - dif[i]) / dif[i]);
		printf("divi error %e\n", altLL);
		altLL = 0.0;
	}

	LL *= -0.5;
	dif = cached_kernel.leftDivision(dif);

	Tuple< Tuple<double> > ders; ders.setSize(cached_kernel.getSize() * scale.getSize());
	Tuple<double> tmp; tmp.setSize(scale.getSize());

	for(m=0;m<ders.getSize();m++) ders[m].setSize(cached_kernel.getSize()).toZero();
	Tuple<uint32_t, 2u> coor, coor2;
	for(k=0,coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
		for(coor2[1]=0;coor2[1]<coor[1];coor2[1]++,k++){ // i and j are sampleIDs
			for(coor[0]=0;coor[0]<scale.getSize();coor[0]++) {coor2[0] = coor[0]; tmp[coor[0]] = timestamp[coor2] - timestamp[coor];}
			//printf("pre\n");tmp.show();
			tmp = scale * tmp; // projection
			//printf("post\n");tmp.show();
			for(m=0;m<scale.getSize();m++){ // n and m are dimensions
				ders[coor[1] * scale.getSize() + m][coor2[1]] -= cached_kernel.data[k] * tmp[m];
				ders[coor2[1] * scale.getSize() + m][coor[1]] += cached_kernel.data[k] * tmp[m];
			}
		}
		k++;
	}

	for(k=0;k<ders.getSize();k++) {
		altLL =  dif[0] * ders[k][0];
		for(m=1;m<cached_kernel.getSize();m++) altLL += dif[m] * ders[k][m];
		altLL *= dif[(k / scale.getSize())];

		ders[k] = cached_kernel.leftDivision(ders[k]);
		altLL = ders[k][(k / scale.getSize())] - altLL;
		fout[k] += 2.0 * altLL;
	}
return LL;}*/



Tuple<double> GaussianProcessMK2::setKernel3(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp) const{
	Tuple<double> dev; dev.setSize(observ.getShape(1));
	cached_kernel.setSize(observ.getShape(1));
	Tuple<double> tmp; tmp.setSize(timestamp.getShape(0));
	Tuple<uint32_t, 2u> coor;
	Tuple<uint32_t, 2u> coor2;
	int k;
	for(k=0,coor[1]=0;coor[1]<observ.getShape(1);coor[1]++){
		coor[0] = 0;
		dev[coor[1]] = observ[coor] - this->getMean();
		for(coor2[1]=0;coor2[1]<coor[1];coor2[1]++){
			for(coor[0]=0;coor[0]<tmp.getSize();coor[0]++) {coor2[0] = coor[0]; tmp[coor[0]] = timestamp[coor2] - timestamp[coor];}
			cached_kernel.data[k++] = fabs(this->getSignal()) * exp( -scale.Xformed_inner_product(tmp));
		}
		coor[0] = 1;
		cached_kernel.data[k++] = fabs(this->getSignal()) + fabs(this->getNoise()) + observ[coor];
	}
	cached_projection = cached_kernel.leftDivision(dev);
return dev;}

double GaussianProcessMK2::getLL(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp){ double fout;
	short positivedefclass;
	Tuple<double> dif; dif = this->setKernel3(observ,timestamp);
	fout = cached_kernel.log_determinant(&positivedefclass);
	if (positivedefclass != 1) return nan("NAN");
return -0.5 * (fout + cached_kernel.Xformed_inner_product_of_inverse(dif));}
double GaussianProcessMK2::gradientAscent(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, unsigned int nbsteps, double *startLL){
	CurvatureSearchScope dacurv;
	dacurv.init();
	Tuple<double> deriv;
	deriv.setSize(3 + scale.totsize());
	Tuple<double> dif; dif.setSize(timestamp.dims[1]);
	Trianglix< Tuple<double> > ders; ders.setSize(dif.getSize());
	Tuple<double> tmp; tmp.setSize(scale.getSize());
	Tuple<uint32_t, 2u> coor;
	Tuple<uint32_t, 2u> coord;
	int k;
	for(k=0,coor[1]=0;coor[1]<dif.getSize();coor[1]++){
		for(coord[1]=0;coord[1]<=coor[1];coord[1]++) {
			ders.data[k].setSize(2 + scale.totsize());
			ders.data[k++][0] = 0.0;
		}
	}
	SparseMatrix<double> tracetrace;tracetrace.setNBcols(10);
/*
	for(coor[1]=0;coor[1]<dif.getSize();coor[1]++){
		coor[0]=0;
		printf("%i[%e", coor[1], timestamp[coor]);
		coor[0]=1;
		printf(", %e]: ", timestamp[coor]);
		coor[0]=0;
		printf("%e +- ",observ[coor]);
		coor[0]=1;
		printf("%e\n", observ[coor]);
	}*/


	uint32_t curstep;
	double tmpLL, tmpLL2;
	short positivedefclass;
	short positivedefclass2;

	for(curstep = 0; curstep<nbsteps;curstep++){
		unsigned int l,m,n;
		dif = this->setKernel3(observ,timestamp);
		tmpLL2 = cached_kernel.Xformed_inner_product_of_inverse(dif);

		tmpLL = cached_kernel.log_determinant(&positivedefclass);
		 // has non-positive eigen values...
		double LL = -0.5 * (tmpLL + tmpLL2);
		scale.log_determinant(&positivedefclass2);

		if ((positivedefclass != 1)||(positivedefclass2 != 1)||(!ExOp::isValid(LL))) {
			if (curstep == 0) {
				this->getSignal() *= 0.125;
				curstep--;
				continue;
			}
			if (dacurv.updateGotNAN(this->mkParamIterator())) break;
			continue;
		}

		//printf("LL[%i] = %e (%e + %e) %e %e %e\n", curstep, -0.5*(LL + tmpLL), -0.5*LL, -0.5*tmpLL, data[0], data[1], data[2]);

		//Tuple<double> old = dif;
		dif = cached_kernel.leftDivision(dif);
		// derivative in Signal
		//Trianglix<double> der, der2; der.setSize(cached_kernel.getSize()); der2.setSize(cached_kernel.getSize());
		deriv[0] = 0.0;
		for(k=0,coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
			for(coord[1]=0;coord[1]<coor[1];coord[1]++,k++){
				for(coord[0]=0;coord[0]<scale.getSize();coord[0]++) {coor[0] = coord[0]; tmp[coord[0]] = timestamp[coord] - timestamp[coor];}
				ders.data[k][1] = -cached_kernel.data[k] / fabs(getSignal());
				//der.data[k] = cached_kernel.data[k] / getSignal();
				//der2.data[k] = 0.0;
				for(l=2,m=0;m<scale.getSize();m++){
					for(n=0;n<m;n++) ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[n] * 2.0;
					ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[m];
				}
			}
			deriv[0] += dif[coor[1]];
			//fout[1] -= noisescale[i] * (dif[i] * dif[i] + 1.0 / cached_kernel.data[k]);
			ders.data[k].toZero();
			ders.data[k][1] = -1.0;
			coor[0] = 1; ders.data[k][0] = -1.0; // observ[coor];
			//der2.data[k] = noisescale[i];
			//der.data[k] = 1.0;
			k++;

		}

		Tuple<double> trofdif = ders.Xformed_inner_product(dif);
		trofdif -= ders.trace_of_division(cached_kernel);
		for(l=0;l< trofdif.getSize();l++) deriv[l+1] = trofdif[l] * -0.5;
		int nberr = 0;
		if (auto dite = deriv.getIterator()) do{
			if (!ExOp::isValid(*dite)) nberr++;
		}while(dite++);

		if ((!ExOp::isValid(LL))||((nberr != 0)&&(LL > dacurv.getLastValue()) )){
			short hehe;
			if ((!ExOp::isValid(tmpLL))&&(ExOp::isValid(tmpLL2))) cached_kernel.log_determinant(&hehe,true);
			printf("Got NAN in deriv %e %e:\n", tmpLL, tmpLL2);
			if (auto ite = this->mkParamIterator()) if (auto dite = deriv.getIterator()) do{
				printf("%e\t%e\n", *ite, *dite);
				ite++;
			}while(dite++);
			dif.show();
			cached_kernel.show();
		}

		//if (curstep < 24){
		if (this->getSignal() < 0.0) deriv[2] = -deriv[2];
		if (this->getNoise() < 0.0) deriv[1] = -deriv[1];

		if (dacurv.updateAscent(LL,this->mkParamIterator(),deriv.getIterator(), &tracetrace , curstep)) break;
		//	printf("%e\t%e\t%e\t%e\t%e\t%e\n", data[0], data[1], data[2], scale.data[0], scale.data[1], scale.data[2]);
		this->checkDomain(dacurv);
		//}else{
		//	if (dacurv.checkDerivative(LL,this->mkParamIterator(),deriv.getIterator())) break;
		//}
		if ((curstep == 0)&&(startLL != NULL)) startLL[0] = LL;
	}
	//tracetrace.show();
	if (curstep == nbsteps) dacurv.wrFinalGuess(this->mkParamIterator());
	if (this->getSignal() < 0.0) this->getSignal() = -this->getSignal();
	if (this->getNoise() < 0.0) this->getNoise() = -this->getNoise();
	this->setKernel3(observ,timestamp);


	double toreturn = dacurv.getLastValue();
	if (!ExOp::isValid(toreturn)) {
		toreturn = 0.0;
		this->show();
		cached_kernel.show();
	}
	double tmpLLval = this->getLL(observ,timestamp);
	if (tmpLLval != toreturn) foreach.printf("LL mismatch (step %i)! %e vs %e\n", curstep, toreturn, this->getLL(observ,timestamp));
return toreturn;}

double GaussianProcessMK2::gradientAscentMK2(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, unsigned int nbsteps){
	AdaptiveLearningScope alsc;

	Tuple<double> deriv;
	deriv.setSize(3 + scale.totsize());
	Tuple<double> dif; dif.setSize(cached_kernel.getSize());
	Trianglix< Tuple<double> > ders; ders.setSize(cached_kernel.getSize());
	Tuple<double> tmp; tmp.setSize(scale.getSize());
	Tuple<uint32_t, 2u> coor;
	Tuple<uint32_t, 2u> coord;
	int k;
	for(k=0,coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
		for(coord[1]=0;coord[1]<=coor[1];coord[1]++) {
			ders.data[k].setSize(2 + scale.totsize());
			ders.data[k++][0] = 0.0;
		}
	}
	SparseMatrix<double> tracetrace;tracetrace.setNBcols(10);

	for(coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
		coor[0]=0;
		printf("%i[%e", coor[1], timestamp[coor]);
		coor[0]=1;
		printf(", %e]: ", timestamp[coor]);
		coor[0]=0;
		printf("%e +- ",observ[coor]);
		coor[0]=1;
		printf("%e\n", observ[coor]);
	}



	double tmpLL;
	if (auto alp= AdativeLearningParameters()) do{
		unsigned int l,m,n;
		dif = this->setKernel3(observ,timestamp);
		alp.cur_value = cached_kernel.Xformed_inner_product_of_inverse(dif);
		tmpLL = cached_kernel.log_determinant();
		printf("LL[%i] = %e (%e + %e) %e %e %e\n", alp.step, -0.5*(alp.cur_value + tmpLL), -0.5*alp.cur_value, -0.5*tmpLL, data[0], data[1], data[2]);
		alp.cur_value += tmpLL;
		alp.cur_value *= -0.5;
		//Tuple<double> old = dif;
		dif = cached_kernel.leftDivision(dif);
		// derivative in Signal
		//Trianglix<double> der, der2; der.setSize(cached_kernel.getSize()); der2.setSize(cached_kernel.getSize());
		deriv[0] = 0.0;
		for(k=0,coor[1]=0;coor[1]<cached_kernel.getSize();coor[1]++){
			for(coord[1]=0;coord[1]<coor[1];coord[1]++,k++){
				for(coord[0]=0;coord[0]<scale.getSize();coord[0]++) {coor[0] = coord[0]; tmp[coord[0]] = timestamp[coord] - timestamp[coor];}
				ders.data[k][1] = -cached_kernel.data[k] / getSignal();
				//der.data[k] = cached_kernel.data[k] / getSignal();
				//der2.data[k] = 0.0;
				for(l=2,m=0;m<scale.getSize();m++){
					for(n=0;n<m;n++) ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[n] * 2.0;
					ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[m];
				}
			}
			deriv[0] += dif[coor[1]];
			//fout[1] -= noisescale[i] * (dif[i] * dif[i] + 1.0 / cached_kernel.data[k]);
			ders.data[k].toZero();
			ders.data[k][1] = -1.0;
			coor[0] = 1; ders.data[k][0] = -1.0; // observ[coor];
			//der2.data[k] = noisescale[i];
			//der.data[k] = 1.0;
			k++;
		}


		Tuple<double> trofdif = ders.Xformed_inner_product(dif);
		trofdif -= ders.trace_of_division(cached_kernel);
		for(l=0;l< trofdif.getSize();l++) deriv[l+1] = trofdif[l] * -0.5;

		if (auto dite = deriv.getIterator()) do{printf("%e\t", *dite);
		}while(dite++);
		printf("\n");
		//if (curstep < 174){
			if (alsc.updateVariable(alp,this->mkParamIterator(),deriv.getIterator())) break;
			this->checkDomain(alsc);
		/*else{

			if (alsc.checkDerivative(LL,this->mkParamIterator(),deriv.getIterator())) break;
		}*/

	}while(alp++);
	//tracetrace.show();
	//if (curstep == nbsteps) dacurv.wrFinalGuess(this->mkParamIterator());
return 0.0;}

DataGrid<double, 3u> GaussianProcessMK2::makeImage(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, int sizex, int sizey, int borderx, int bordery, Tuple<double> *extra) const{DataGrid<double, 3u> fout;
	Tuple<uint32_t, 3u> coor; coor[0] =2; coor[1] =sizex; coor[2] =sizey; fout.setSizes(coor);
	class SharedScope{ public:
		const GaussianProcessMK2& gp;
		const RemoteMemory<const double,2u> &observ;
		const RemoteMemory<const double,2u> &timestamp;
		Tuple<double> projection;
		//double signalandnoise;
		Tuple<Tuple<double,2u>, 0u> datarange;
		SharedScope(const GaussianProcessMK2& _gp, const RemoteMemory<const double,2u> &_observ, const RemoteMemory<const double,2u> &_timestamp):gp(_gp), observ(_observ), timestamp(_timestamp){}
	}shrd(*this,observ, timestamp);

	class Task{
	public:
		double range[4];
		Task(){range[0] = 1.0; range[1]= 0.0;}

		void operator()(AtomicIterator<decltype(DataGrid<double, 3u>().getPixelIterator())> &pixels, const SharedScope& shrd){
			Tuple<double> dif; dif.setSize(shrd.gp.scale.getSize());
			unsigned int i,l;
			Tuple<double> position; position.setSize(shrd.gp.scale.getSize());
			Tuple<uint32_t, 2u> coor; coor[0] =0;
			Tuple<double> statedif; statedif.setSize(shrd.gp.scale.getSize());
			Tuple<double> localkernel; localkernel.setSize(shrd.gp.cached_kernel.getSize());

			//foreach.print("%X is address of doubles %i", (intptr_t)localkernel(), dif.getSize());

			if (auto pite = pixels()) do{
				// x ~ N( v K^-1 f,  v K^-1 v + noise)
				if (auto tite = shrd.timestamp.mkIterator()){
					for(l=0;l<localkernel.getSize();l++) {
						for(i=0;i<dif.getSize();i++, tite++) dif[i] = shrd.datarange[i][1] * pite()[i] + shrd.datarange[i][0] - (*tite);
						//foreach.print("%e vs %e haha\n", dif[0]*dif[0] + dif[1]*dif[1], shrd.gp.scale.Xformed_inner_product(dif));
						localkernel[l] = shrd.gp.getSignal() * exp( -shrd.gp.scale.Xformed_inner_product(dif));
					}
				}

				//localkernel.toZero();
				//localkernel[0] = shrd.datarange[0][1] * pite()[0] + shrd.datarange[0][0];
				//localkernel[1] = shrd.datarange[1][1] * pite()[1] + shrd.datarange[1][0];

				(*pite)[0] = localkernel.mkInnerProd(shrd.projection); // expected mean
				(*pite)[1] = shrd.gp.cached_kernel.Xformed_inner_product_of_inverse(localkernel); // explained variance

				//(*pite)[0] = shrd.gp.cached_kernel.Xformed_inner_product_of_inverse(localkernel); //shrd.datarange[0][1] * pite()[0] + shrd.datarange[0][0];
				//(*pite)[1] = localkernel[0]; //shrd.datarange[1][1] * pite()[1] + shrd.datarange[1][0];

				if (range[0] >range[1] ) {
					range[0] = range[1] = (*pite)[0];
					range[2] = range[3] = (*pite)[1];
				}else{
					if ((*pite)[0] < range[0]) range[0] = (*pite)[0];
					else if ((*pite)[0] > range[1]) range[1] = (*pite)[0];
					if ((*pite)[1] < range[2]) range[2] = (*pite)[1];
					else if ((*pite)[1] > range[3]) range[3] = (*pite)[1];
				}
			}while(pite++);
			//foreach.printf("done\n");
		}
		void operator<<=(Task& other){
			if (range[0] >range[1]) {
				range[0] = other.range[0];
				range[1] = other.range[1];
				range[2] = other.range[2];
				range[3] = other.range[3];
			}else if (other.range[0] <= other.range[1]) {
				if (other.range[0]< range[0]) range[0] = other.range[0];
				else if (other.range[1] > range[1]) range[1] = other.range[1];
				if (other.range[2]< range[2]) range[2] = other.range[2];
				else if (other.range[3] > range[3]) range[3] = other.range[3];
			}
		}
		Tuple<double, 4u> operator()(const SharedScope& shrd){ Tuple<double, 4u> fout;
			fout[0] = range[0];
			fout[1] = range[1];
			fout[2] = range[2];
			fout[3] = range[3];
		return fout;}
	};

	if (observ.getShape()[0] != 2) {printf("Observation shape[0] is expected to be 2, obs and stddev interleaved"); return fout;}
	if (timestamp.getShape()[0] != scale.getSize()) {printf("state shape[0] is expected to match state size of GP (which is %i)", scale.getSize()); return fout;}


	// Step 1, solve K y = projection, and get coordinate bounds

	shrd.projection = this->cached_kernel.leftDivision(this->setKernel3(observ,timestamp));


	shrd.datarange.setSize(timestamp.getShape(0)); // should be 2 ...
	if (auto ite = timestamp.mkIterator()) {
		shrd.datarange[0][0] = shrd.datarange[0][1] = *ite;
		for(unsigned int i=1; i < shrd.datarange.getSize();i++){
			ite++; shrd.datarange[i][0] = shrd.datarange[i][1] = *ite;
		}
		while(ite++){
			if (*ite < shrd.datarange[0][0]) shrd.datarange[0][0] = *ite;
			else if (*ite > shrd.datarange[0][1]) shrd.datarange[0][1] = *ite;
			for(unsigned int i=1; i < shrd.datarange.getSize();i++){
				ite++;
				if (*ite < shrd.datarange[i][0]) shrd.datarange[i][0] = *ite;
				else if (*ite > shrd.datarange[i][1]) shrd.datarange[i][1] = *ite;
			}
		};
	}
	//double harmonic =0.0;
	//Tuple<uint32_t, 2u> coor2; coor2[0] =1;
	//for(coor2[1] =0; coor2[1]< observ.getShape(1);coor2[1]++) harmonic += 1.0 / observ[coor2];
	//printf("harmonic = %e\n", harmonic / observ.getShape(1));
	//shrd.signalandnoise = this->getSignal() + (this->getNoise() * observ.getShape(1)) / harmonic;
	shrd.datarange[0][1] = (shrd.datarange[0][1] - shrd.datarange[0][0]) / (sizex - 1 - borderx);
	shrd.datarange[1][1] = (shrd.datarange[1][1] - shrd.datarange[1][0]) / (sizey - 1 - bordery);
	shrd.datarange[0][0] -= (shrd.datarange[0][1] * (borderx >>1));
	shrd.datarange[1][0] -= (shrd.datarange[1][1] * (bordery >>1));

	//shrd.gp.setKernel3();
	Tuple<double, 4u> res = foreach[fout.getPixelIterator()].run<Task>(shrd);
    if (extra != NULL) {
		extra->setSize(8);
		(*extra)[0] = shrd.datarange[0][0];
		(*extra)[1] = shrd.datarange[0][1];
		(*extra)[2] = shrd.datarange[1][0];
		(*extra)[3] = shrd.datarange[1][1];
		(*extra)[4] = res[0];
		(*extra)[5] = res[1];
		(*extra)[6] = res[2];
		(*extra)[7] = res[3];
    }
return fout;}

GaussianProcessArray& GaussianProcessArray::initUsingPCA(RemoteMemory<const double,2u> observ, uint32_t nbPCs){
	printf("init with %iPCs, from %i x %i matrix\n", nbPCs, observ.dims[0], observ.dims[1]);
	hstate.setSizes(nbPCs ,observ.dims[0]);
	RemoteMemory<const double,2u> obsonly = observ.toSwapDims(0,1);
	TMatrix<double> upper; (upper = obsonly).toUpperTriangular();

	Tuple<uint32_t , 2u> coor, coor2;
	double fact;
	coor[0] = 0;
	coor[1] = 0;
	fact = upper(coor);
	for(coor[1]++;coor[1] < upper.sizes[1];coor[1]++) fact += upper(coor);
	fact /= upper.sizes[1];
	for(coor[1]=0;coor[1] < upper.sizes[1];coor[1]++) upper(coor) -= fact;
	for(coor[0]++;coor[0]<upper.sizes[1];coor[0]++){
		coor[1] = coor[0]; fact = upper(coor);
		for(coor[1]++;coor[1] < upper.sizes[1];coor[1]++) fact += upper(coor);
		fact /= upper.sizes[1];
		for(coor[1] = 0; coor[1] < upper.sizes[1];coor[1]++) upper(coor) -= fact;
	}

	Tuple<double> eigenval;
	TMatrix<double> eigenvec = upper.mkOuterProduct().getEigenVectors(eigenval);
	if (auto hite = hstate.getIterator()) do{*hite = eigenvec(hite());}while(hite++);
return *this;}

double GaussianProcessArray::cacheObservations(RemoteMemory<const double,3u> observ){
	RemoteMemory<double,2u> mataccess = hstate;
class GlobalScope{ public:
	RemoteMemory<const double,3u>& observ;
	RemoteMemory<double,2u>& timestamp;
	mutable bool gotNAN;
	bool doverbose;

	GlobalScope(RemoteMemory<const double,3u>& _observ, RemoteMemory<double,2u>& _timestamp): observ(_observ), timestamp(_timestamp){}
}shrd(observ, mataccess);
class TaskCacheProj{ public:
	double LL;
	TaskCacheProj(): LL(0.0){}
	void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
		double tmp;
		if (auto jobite = cgp()) do{
			 tmp = jobite->getLL(shrd.observ[jobite()],shrd.timestamp);
			 if (ExOp::isValid(tmp)) LL += tmp;
			 else foreach.print("Nan for feature %i\n", jobite());
		}while(jobite++);
	}
	void operator<<=(const TaskCacheProj& other){LL += other.LL;}
	double operator()(const GlobalScope& shrd){return LL;}
};
return foreach[gps.mkIterator()].run<TaskCacheProj>(shrd);}

double GaussianProcessArray::learnThis(RemoteMemory<const double,3u> observ, uint32_t nbsteps){
	if (nbsteps == 0) return this->cacheObservations(observ);
	if (observ.dims[0] != 2) {printf("expects tuples of length 2, observation and std err\n"); return 0.0;}
	if (gps.getSize() != observ.dims[2]) {
		if (gps.getSize() != 0){
			printf("warning: number of observation differs, expected %i but got %i.\n", gps.getSize(), observ.dims[1]);
			this->setSize(observ.dims[1],gps[0].scale.getSize());
		}else{
			printf("warning: setting default size 2.\n");
			this->setSize(observ.dims[1],2);
		}
	}
	if (hstate.sizes[1] != observ.dims[1]) this->initUsingPCA(observ, gps[0].scale.getSize());


	RemoteMemory<double,2u> mataccess = hstate;
	class GlobalScope{ public:
		RemoteMemory<const double,3u>& observ;
		RemoteMemory<double,2u>& timestamp;
		mutable bool gotNAN;
		bool doverbose;

		GlobalScope(RemoteMemory<const double,3u>& _observ, RemoteMemory<double,2u>& _timestamp): observ(_observ), timestamp(_timestamp){}
	}shrd(observ, mataccess);

	class TaskGPparam{ public:
		double LL;
		TaskGPparam(): LL(0.0){}
		void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
			double start,tmp;
			if (auto jobite = cgp()) do{
				 tmp = jobite->gradientAscent(shrd.observ[jobite()],shrd.timestamp,1000, &start);
				 if (ExOp::isValid(tmp)) LL += tmp;
				 else foreach.printf("Delta[%i]: %e (%e)\n", jobite(), tmp - start, tmp);
			}while(jobite++);
		}
		void operator<<=(const TaskGPparam& other){LL += other.LL;}
		double operator()(const GlobalScope& shrd){return LL;}
	};

	class TaskStateParam{ public:
		double LL;
		Tuple<double> local_deriv;
		TaskStateParam():LL(0.0){}
		void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
			LL = 0.0;
			local_deriv.setSize(shrd.timestamp.dims[0] * shrd.timestamp.dims[1]).toZero();
			double tmp;
			if (auto jobite = cgp()) do{
				tmp += jobite->addLLTimeDerivative(local_deriv,shrd.observ[jobite()],shrd.timestamp, (shrd.doverbose)&&(jobite() == 0));
				if (!ExOp::isValid(tmp)) {
					shrd.gotNAN = true;
					foreach.printf("Got nan for %i\n", jobite());
					foreach.show(*jobite);
					foreach.show(jobite->cached_kernel);
				}else LL += tmp;
			}while(jobite++);
			//foreach.printf("%e",LL);
		}
		void operator<<=(const TaskStateParam& other){LL += other.LL; local_deriv += other.local_deriv;}
		KeyElem<double, Tuple<double>>  operator()(const GlobalScope& shrd){KeyElem<double, Tuple<double>> fout; fout.k = LL; fout.d.toMemmove(local_deriv); return fout;}
	};


	CurvatureSearchScope dacsscp;
	KeyElem<double, Tuple<double>> res;
	int innerstep;
	SparseMatrix<double> logout; logout.setNBcols(10);
	shrd.doverbose = false;
	for(uint32_t j=0;j<nbsteps;j++){
	    double dall = foreach[gps.mkIterator()].run<TaskGPparam>(shrd);
	    printf("daLL after %e\n",dall);
		dacsscp.init();
		for(innerstep=0;innerstep < 20; innerstep++){
			shrd.gotNAN = false;
			res = foreach[gps.mkIterator()].run<TaskStateParam>(shrd);

			/*if (auto ite = timestamp.mkIterator()) if (auto dite = res.d.mkIterator()) do{
				printf("%i:%e\t%e\n", k++, *ite, *dite);
				ite++;
			}while(dite++);*/
			if (shrd.gotNAN) res.k = sqrt(-1.0);

			if (innerstep < 20){
				printf("LL[%i] = %e\n", innerstep, res.k);
				if (dacsscp.updateAscent(res.k, hstate.mkIterator(), res.d.mkIterator(),&logout,innerstep)) break;
			}else {
				if (innerstep == 20) {dacsscp.init(); shrd.doverbose = true;}
				dacsscp.checkDerivative(res.k, hstate.mkIterator(), res.d.mkIterator());
			}
		}
		if (innerstep == 1000) dacsscp.wrFinalGuess(hstate.mkIterator());
		logout.show();
	}
return 0.0;}
double GaussianProcessArray::learnGPParams(RemoteMemory<const double,3u> observ){
	if (hstate.data.getSize() == 0) {printf("expects an initialized hidden state\n"); return 0.0;}
	if (hstate.sizes[1] != observ.dims[1]) {printf("mitmatch for number of observation %i, hidden states implies %i\n", observ.dims[1], hstate.sizes[1]); return 0.0;}
	if (observ.dims[0] != 2) {printf("expects tuples of length 2, observation and std err\n"); return 0.0;}
	gps.setSize(observ.dims[2]);
	for(unsigned int i=0;i<gps.getSize();i++) gps[i].init();
	RemoteMemory<double,2u> mataccess = hstate;
	class GlobalScope{ public:
		RemoteMemory<const double,3u>& observ;
		RemoteMemory<double,2u>& timestamp;
		mutable bool gotNAN;
		bool doverbose;
		GlobalScope(RemoteMemory<const double,3u>& _observ, RemoteMemory<double,2u>& _timestamp): observ(_observ), timestamp(_timestamp){}
	}shrd(observ, mataccess);
	class TaskGPparam{ public:
		double LL;
		TaskGPparam(): LL(0.0){}
		void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
			double start,tmp;
			if (auto jobite = cgp()) do{
				 tmp = jobite->gradientAscent(shrd.observ[jobite()],shrd.timestamp,1000, &start);
				 if (ExOp::isValid(tmp)) LL += tmp;
				 else foreach.printf("Delta[%i]: %e (%e)\n", jobite(), tmp - start, tmp);
				 if ((jobite->cached_projection.getSize() != jobite->cached_kernel.getSize())||(!ExOp::isValid(jobite->cached_projection))){
					foreach.print("Learn failed for %i\n",jobite());
					foreach.show(jobite->cached_projection);
				 }
			}while(jobite++);
		}
		void operator<<=(const TaskGPparam& other){LL += other.LL;}
		double operator()(const GlobalScope& shrd){return LL;}
	};
	CurvatureSearchScope dacsscp;
	KeyElem<double, Tuple<double>> res;
	int innerstep;
	SparseMatrix<double> logout; logout.setNBcols(10);
	shrd.doverbose = false;
	double dall = foreach[gps.mkIterator()].run<TaskGPparam>(shrd);
return 0.0;}
TMatrix<double> GaussianProcessArray::learnPCAreduction(RemoteMemory<const double,4u> observ, uint32_t nb_latent, uint32_t nbsteps, uint32_t nbpc_per_ct){
	if (observ.dims[0] != 2) {printf("expects tuples of length 2, observation and std err\n"); return TMatrix<double>();}
	if (nbpc_per_ct == 0) nbpc_per_ct = observ.dims[1] -1;
	else if (nbpc_per_ct >= observ.dims[1]) nbpc_per_ct = observ.dims[1] - 1;
	RemoteMemory<double,2u> mataccess = hstate;
	RemoteMemory<const double,3u> thredaccess = observ.vectorizeOntoDimension(3);
class GlobalScope{ public:
	RemoteMemory<const double,4u>& observ;
	RemoteMemory<const double,3u>& observ3;
	RemoteMemory<double,2u>& timestamp;
	mutable bool gotNAN;
	bool doverbose;
	mutable Tuple<double> derivative;
	DataGrid<double, 4u> reduced; // 2 x (Donor) x (nbpc_per_ct) x (Celltypers)
	TMatrix<double> fout;
	uint32_t gpderivsize;
	GlobalScope(RemoteMemory<const double,4u>& _observ, RemoteMemory<double,2u>& _timestamp,RemoteMemory<const double,3u>& _observ3, uint32_t nbpc_per_ct): observ(_observ), timestamp(_timestamp),observ3(_observ3){
		gpderivsize = 3 + ((timestamp.dims[0] * (timestamp.dims[0] +1)) >> 1);
		derivative.setSize(timestamp.dims[0] * timestamp.dims[1] + gpderivsize * observ.dims[2]);
		fout.setSizes(nbpc_per_ct * observ.dims[3], observ.dims[2]);
		Tuple<uint32_t, 4u> dim; dim[0] = 2; dim[1] = observ.dims[1]; dim[2] = nbpc_per_ct; dim[3] = observ.dims[3];
		reduced.setSizes(dim);

		// to see weird things
		reduced.toZero();
		fout.toZero();
	}
	uint32_t getNBcelltypes() const{return observ.dims[3];}
}shrd(observ, mataccess, thredaccess, nbpc_per_ct);

class PCAonInput{ public:
	PCAonInput(){}
	void operator<<=(const PCAonInput& other){}
	double operator()(GlobalScope& shrd){return 0.0;}
	void operator()(AtomicIterator<RangeIterator>& cgp, GlobalScope& shrd){
		if (auto jobite = cgp()) do{
			RemoteMemory<const double,2u> obsonly = shrd.observ.selectSlice(jobite(),3).selectSlice(0,0).toSwapDims(0,1);
			TMatrix<double> ortho;
			Tuple<double> eigen;
			TMatrix<double> topQ; topQ.toSingularDecomposition(obsonly, eigen ,ortho);
			// now topQ * obsonly is a n x n matrix
			Tuple<uint32_t,2> coor2, coor2alt;
			Tuple<uint32_t,4> coor4;
			coor4[3] = jobite();
			RemoteMemory<const double,2u> unsonly = shrd.observ.selectSlice(jobite(),3).selectSlice(1,0);

			TMatrix<double> daprojchk = topQ * shrd.observ.selectSlice(jobite(),3).selectSlice(0,0).toSwapDims(0,1);

			foreach.show(daprojchk);
			if (auto hite = ortho.getIterator()) do{ // tr(ortho) * diag(eigen) is the projection
				if (hite()[0] < shrd.reduced.dims[2]){
					coor4[1] = hite()[1]; // all donors
					coor4[2] = hite()[0]; // topQ vectors

					//shrd.reduced(coor4) = *hite * eigen[coor4[2]];
					coor4[0] = 0;
					coor2[0] = hite()[0];
					coor2[1] = hite()[1];
					shrd.reduced(coor4) = daprojchk(coor2);

					coor4[0] = 1;
					coor2[1] = 0; coor2alt[1] = 0; coor2[0] = hite()[1]; coor2alt[0] = hite()[0];

					shrd.reduced(coor4) = unsonly[coor2] * pow(topQ(coor2alt), 2.0);
					for(coor2[1]++; coor2[1] < unsonly.dims[1];coor2[1]++) {coor2alt[1] = coor2[1]; shrd.reduced(coor4) += unsonly[coor2] * pow(topQ(coor2alt), 2.0);}
				}
			}while(hite++);


			// define projection...!!!
			if (auto hite = topQ.getIterator()) do{
				coor2[1] = hite()[1];
				if (hite()[0] < shrd.reduced.dims[2]){
					coor2[0] = hite()[0] + jobite() * shrd.reduced.dims[2];
					shrd.fout(coor2) = *hite;
				}
			}while(hite++);
		}while(jobite++);
	}
};

	class TaskAllParam{ public:
		double LL;
		Tuple<double> local_deriv;
		TaskAllParam():LL(0.0){}
		void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
			LL = 0.0;
			unsigned int k,l,m,n;
			Tuple< Tuple<double> > ders; ders.setSize(shrd.timestamp.dims[0] * shrd.timestamp.dims[1]);
			Tuple<uint32_t, 2u> coor, coor2;
			Trianglix< Tuple<double> > gpders; gpders.setSize(shrd.timestamp.dims[1]);

			for(k=0,m=0;m<shrd.timestamp.dims[1];m++,k++){
				for(n=0;n<m;n++,k++) {gpders.data[k].setSize(shrd.gpderivsize-1); gpders.data[k][0] = 0.0;}
				gpders.data[k].setSize(shrd.gpderivsize-1).toZero();
				gpders.data[k][0] = -1.0;
				gpders.data[k][1] = -1.0;
			}

			Tuple<double> tmp; tmp.setSize(shrd.timestamp.dims[0]);
			local_deriv.setSize(shrd.timestamp.dims[0] * shrd.timestamp.dims[1]).toZero();
			for(m=0;m<ders.getSize();m++) ders[m].setSize(shrd.timestamp.dims[1]);
			if (auto jobite = cgp()) do{
				double *gpderiv = shrd.derivative() + jobite() * shrd.gpderivsize;
				Tuple<double> dif = jobite->setKernel3(shrd.observ3[jobite()], shrd.timestamp);
				short posdefclass;
				double altLL = jobite->cached_kernel.Xformed_inner_product_of_inverse(dif);
				LL += jobite->cached_kernel.log_determinant(&posdefclass);
				if (posdefclass != 1) {shrd.gotNAN = true; break;}
				altLL = 0.0;

				if (verbose){
					Tuple<double> andback = jobite->cached_kernel.leftDivision(dif);
					for(unsigned int i =0;i<dif.getSize();i++) altLL += andback[i] * dif[i];
					printf("LL %e vs %e  %e error\n", LL, altLL, (LL - altLL) /LL);
					andback = jobite->cached_kernel * andback;
					altLL = 0.0;
					for(unsigned int i =0;i<dif.getSize();i++) altLL += fabs((andback[i] - dif[i]) / dif[i]);
					printf("divi error %e\n", altLL);
					altLL = 0.0;
				}

				LL *= -0.5;
				dif = jobite->cached_kernel.leftDivision(dif);

				for(m=0;m<ders.getSize();m++) ders[m].toZero();

				gpderiv[0] = 0.0; gpders.toZero();
				for(k=0,coor[1]=0;coor[1]<dif.getSize();coor[1]++,k++){
					gpderiv[0] += dif[coor[1]];
					for(coor2[1]=0;coor2[1]<coor[1];coor2[1]++,k++){ // i and j are sampleIDs
						for(coor[0]=0;coor[0]<tmp.getSize();coor[0]++) {coor2[0] = coor[0]; tmp[coor[0]] = shrd.timestamp[coor2] - shrd.timestamp[coor];}
						tmp = jobite->scale * tmp; // projection
						gpders.data[k][1] = -jobite->cached_kernel.data[k] / jobite->getSignal();
						for(l=2,m=0;m<tmp.getSize();m++){ // n and m are dimensions
							ders[coor[1] * shrd.timestamp.dims[0]+ m][coor2[1]] -= jobite->cached_kernel.data[k] * tmp[m];
							ders[coor2[1] * shrd.timestamp.dims[0] + m][coor[1]] += jobite->cached_kernel.data[k] * tmp[m];
							for(n=0;n<m;n++) gpders.data[k][l++] = jobite->cached_kernel.data[k] * tmp[m] * tmp[n] * 2.0;
							gpders.data[k][l++] = jobite->cached_kernel.data[k] * tmp[m] * tmp[m];
						}
					}
				}

				Tuple<double> trofdif = gpders.Xformed_inner_product(dif);
				trofdif -= gpders.trace_of_division(jobite->cached_kernel);
				for(l=0;l< trofdif.getSize();l++) gpderiv[l+1] = trofdif[l] * -0.5;

				for(k=0;k<ders.getSize();k++) {
					altLL =  dif[0] * ders[k][0];
					for(m=1;m<dif.getSize();m++) altLL += dif[m] * ders[k][m];
					altLL *= dif[(k / tmp.getSize())];
					ders[k] = jobite->cached_kernel.leftDivision(ders[k]);
					altLL = ders[k][(k / tmp.getSize())] - altLL;
					local_deriv[k] += 2.0 * altLL;
				}
			}while(jobite++);

		}
		void operator<<=(const TaskAllParam& other){LL += other.LL; local_deriv += other.local_deriv;}
		//double  operator()(){return LL;}
		double  operator()(const GlobalScope& shrd){KeyElem<double, Tuple<double>> fout; fout.k = LL; fout.d.toMemmove(local_deriv); return LL;}
	};
class TaskGPparam{ public:
	double LL;
	TaskGPparam(): LL(0.0){}
	void operator()(AtomicIterator<decltype(declval<Tuple<GaussianProcessMK2> >().mkIterator())>& cgp, const GlobalScope& shrd){
		double start,tmp;
		if (auto jobite = cgp()) do{
			 foreach.show(*jobite);
			 foreach.show((*shrd.reduced).vectorizeOntoDimension(2).selectSlice(jobite(),2).getShape());
			 tmp = jobite->gradientAscent((*shrd.reduced).vectorizeOntoDimension(2).selectSlice(jobite(),2),shrd.timestamp,1000, &start);
			 if (ExOp::isValid(tmp)) LL += tmp;
		}while(jobite++);
	}
	void operator<<=(const TaskGPparam& other){LL += other.LL;}
	double operator()(const GlobalScope& shrd){return LL;}
};
	// STEP 1: Reduce Dimensions in all celltypes

	printf("step 1:\n"); fflush(stdout);
	foreach[RangeIterator(observ.dims[3])].run<PCAonInput>(shrd);

	printf("DONE! here is the projections:\n");
	(*shrd.reduced).selectSlice(0,0).show();
	printf("DONE! here is the uncertainty projections:\n");
	(*shrd.reduced).selectSlice(1,0).show();

	printf("DONE! plain projection:\n");

	// STEP 2: Init hidden state with projections
	printf("step 2:\n"); fflush(stdout);
	this->initUsingPCA((*shrd.reduced).selectSlice(0,0).vectorizeOntoDimension(1), nb_latent);
	printf("step 2 done:\n"); fflush(stdout);
	hstate.show();

	// STEP 3: Learn hidden state

	printf("step 3:\n"); fflush(stdout);
	AdaptiveLearningScope dacsscp;
	AdativeLearningParameters daparam;
	this->setSize(shrd.reduced.dims[2]*shrd.reduced.dims[3]);
	for(int i=0;i<gps.getSize();i++) gps[i].init();
	printf("allocated the gps\n");

	for(uint32_t j=0;j<nbsteps;j++){ // well well lets ignore this for now... D:
	    break; // >  _ <
	    double dall = foreach[gps.mkIterator()].run<TaskGPparam>(shrd);

		KeyElem<double, Tuple<double>> res;
		int innerstep;
		SparseMatrix<double> logout; logout.setNBcols(10);
		shrd.doverbose = false;
		uint32_t nbfix;
		if (auto loopite = AdativeLearningParameters()) do{
			loopite.cur_value = foreach[gps.mkIterator()].run<TaskAllParam>(shrd);
			printf("LL[%i] = %e\n", loopite(), loopite.cur_value);
			if (dacsscp.updateVariable(loopite, this->mkIterator(), shrd.derivative.mkIterator())) break;

			nbfix = 0;
			for(unsigned int i=0;i<shrd.timestamp.dims[1];i++){
				if (dacsscp.makeGuessStrictlyPositive(gps[i].data[1], i * shrd.gpderivsize + 1)) {nbfix++; gps[i].data.show();}
				if (dacsscp.makeGuessStrictlyPositive(gps[i].data[2], i * shrd.gpderivsize + 2)) {nbfix++; gps[i].data.show();}

			}
			printf("nbfix = %i\n", nbfix);
		}while(loopite++);
	}
return shrd.fout;}


Tuple< TMatrix<double>, 3u > GaussianProcessArray::localPoissonRegression(RemoteMemory<const uint32_t,2u> observ) const{
class GlobalScope{ public:
	RemoteMemory<const uint32_t,2u>& observ;
	const GaussianProcessArray& gpa;
	uint32_t nbcelltypes;
	Tuple< TMatrix<double>, 3u > fout;
	TMatrix<double> gpcoor;
	Tuple<double> mincoor;
	Tuple<double> offcoor;
	uint32_t gridsize;
	WeightElem<double, 2u> LLspread_stats;
	Tuple<double> LL_out;
	Tuple<uint32_t > Step_out;
	Tuple<uint32_t> Redo_out;
	int nbgenes;
	GlobalScope(const GaussianProcessArray& _gpa, RemoteMemory<const uint32_t,2u> &_observ) : gpa(_gpa), observ(_observ){
		nbgenes = observ.dims[0];
		if (gpa.gps.getSize() == 0) {printf("error, gp not initialized\n");}
		else{
			nbcelltypes = gpa.gps.getSize() / nbgenes;
			fout[0].setSizes(gpa.gps[0].scale.getSize() +2,observ.dims[1]);
			printf("%i vs %i implies %i exposures to infer\n", nbgenes, gpa.gps.getSize(), gpa.gps.getSize() / nbgenes);
			fout[1].setSizes(gpa.gps.getSize() / nbgenes,observ.dims[1]);
			fout[2].setSizes(nbgenes,observ.dims[1]);
			LL_out.setSize(observ.dims[1]);
			Step_out.setSize(observ.dims[1]);
			Redo_out.setSize(observ.dims[1]);
		}
	}
}shrd(*this, observ);
class TaskDelinate{ public:
	double LL, LL_std;
	Tuple<double> local_deriv;
	TaskDelinate():LL(0.0){}
	void operator()(AtomicIterator<decltype(declval<RemoteMemory<const uint32_t,2u> >().mkVectorIterator())>& cgp, GlobalScope& shrd){
		TMatrix<double> LLmatrix; LLmatrix.setSizes(8,8);
		Tuple<double> gpcoor; gpcoor.setSize(shrd.mincoor.getSize()+1);
		Tuple<double> ctmix; ctmix.setSize(shrd.nbcelltypes);
		Tuple<double> guess, deriv; guess.setSize(3 + shrd.nbcelltypes); deriv.setSize(guess.getSize());

		Tuple< uint32_t, 2u> coor;
		KeyElem<double, Tuple<double> > bestguess;
		Tuple<uint32_t, 2u> bestwhich;
		WeightElem<double,2u> LL_stat;
		CurvatureSearchScope dacurv;
		double curLL, LLconstant;
		const int maxstep = 1000;
		uint32_t ascloop;
		double mean, dmean, tmpLL;
		TMatrix<double> deltastate; deltastate.setSizes(shrd.gpa.hstate.sizes[0], shrd.gpa.hstate.sizes[1]);
		Tuple<double> computation_buffer, comp_deriv; computation_buffer.setSize(shrd.gpa.hstate.sizes[1]); comp_deriv.setSize(3);
		Tuple<double> ret_mean; ret_mean.setSize(shrd.nbcelltypes);
		myHashmap<uint32_t> errors;
		if (auto jobite = cgp()) do{

			for(int restart=0;restart<30;restart++){
				guess.toOne();
				guess[0] = shrd.mincoor[0] + (shrd.offcoor[0] * rand()) / RAND_MAX;
				guess[1] = shrd.mincoor[1] + (shrd.offcoor[1] * rand()) / RAND_MAX;

				dacurv.init();
				for(ascloop=0;ascloop<maxstep;ascloop++){
					curLL =0; deriv.toZero();errors.toMemfree();
					if (auto ite = deltastate.getIterator()) do{ *ite = guess[ite()[0]] - shrd.gpa.hstate(ite());} while(ite++);
					for(int g = 0 ; g < shrd.nbgenes; g++){
						mean = 0.0; comp_deriv.toZero();
						for(int c =0; c< shrd.nbcelltypes+3;c++) if (!ExOp::isValid(guess[c])) guess[c] = 1.0;
						for(int c =0; c< shrd.nbcelltypes;c++){ //
							ret_mean[c] = shrd.gpa[c * shrd.nbgenes + g].getMeanAndDeriv_p1(guess[2], exp(guess[c+3]),comp_deriv, deltastate, computation_buffer());
							mean += ret_mean[c];
							/*if ((jobite()[0] == 0)&&(!ExOp::isValid(ret_mean[c]))){
								printf("error for gene%i celltype %i, got %e so we have dagp:\n", g, c, ret_mean[c]);
								shrd.gpa[c * shrd.nbgenes + g].show();
								printf("its scale is:\n");
								shrd.gpa[c * shrd.nbgenes + g].scale.show();
								printf("and projection:\n");
								shrd.gpa[c * shrd.nbgenes + g].cached_projection.show();
								for(int dadada = 0 ; dadada < deltastate.sizes[1]; dadada++) {
									printf("%e -> %e and %e\n", shrd.gpa[c * shrd.nbgenes + g].scale.Xformed_inner_product((*deltastate).selectSlice(dadada,1)), exp(-shrd.gpa[c * shrd.nbgenes + g].scale.Xformed_inner_product((*deltastate).selectSlice(dadada,1))), exp(-guess[2] * guess[2]));
									(*deltastate).selectSlice(dadada,1).show();
									(shrd.gpa[c * shrd.nbgenes + g].scale * (*deltastate).selectSlice(dadada,1)).show();
								}
							}*/
						}
						if ((jobite()[0] == 0)&&(!ExOp::isValid(ret_mean))) {
							//printf("error for gene Count[%i] = %i: \n", g, (*jobite)[g]);
							//ret_mean.show();
							/*if (auto ite = deltastate.getIterator()) do{
								*ite = guess[ite()[0]] - shrd.gpa.hstate(ite());
								printf("%i -> %e; vs %e\n",ite()[0],guess[ite()[0]], shrd.gpa.hstate(ite()) );
								ite().show();
							} while(ite++);
							printf("really?\n");*/
						}
						tmpLL = log(mean) * (*jobite)[g] - mean;
						if (ExOp::isValid(tmpLL)) curLL += tmpLL;
						else if (ascloop == 0) foreach.print("gene %i gives NAs\n", g);
						dmean = ((double)(*jobite)[g]) / mean - 1.0; // dLL / mean[g]
						//if (jobite()[0] == 0) printf("%i: %i vs %e %e\n", (*jobite)[g], mean, dmean);
						if (ExOp::isValid(dmean)) {
							if (ExOp::isValid(comp_deriv))	for(uint32_t cc=0;cc<3;cc++) deriv[cc] += comp_deriv[cc] * dmean;
							if (ExOp::isValid(ret_mean)) for(uint32_t c =0; c< shrd.nbcelltypes;c++) deriv[c+3] += ret_mean[c] * dmean;
						}
						if ((!ExOp::isValid(comp_deriv))||(!ExOp::isValid(ret_mean))||(!ExOp::isValid(dmean))) errors[g] = 1.0;
					}
					/*if (jobite()[0] == 0){
						//errors.show();
						printf("state at %i:\n",ascloop);
						if (auto ite = guess.mkIterator()) do{
							printf("%i:%e,%e\n", ite(), *ite, deriv[ite()]);
						}while(ite++);

						if (dacurv.checkDerivative(curLL,guess.mkIterator(),deriv.mkIterator())) break;
					}else{*/
						if (dacurv.updateAscent(curLL,guess.mkIterator(),deriv.mkIterator())) break;

					//}
				}
				curLL = dacurv.getLastValue();
				//LL_stat += WeightElem<double,2u>(LL);
				if ((restart == 0)||(bestguess.k < curLL)){
					if (ascloop == maxstep) dacurv.wrFinalGuess(guess.mkIterator());
					else ascloop = ascloop - dacurv.getNbStalledSteps();
					bestguess.k = curLL;
					bestguess.d = guess;
					bestwhich[0] = ascloop;
					bestwhich[1] = restart;
				}
			}


			// add constant term to llikelihood
			for(int g = 0 ; g < shrd.nbgenes; g++) bestguess.k -= lngamma(1.0 + (*jobite)[g]);

			coor[1] = jobite()[0];
			for(coor[0]=0;coor[0]< deltastate.sizes[0]+2;coor[0]++)	shrd.fout[0](coor) = (coor[0] < 3) ? bestguess.d[coor[0]] : bestguess.k;
			for(coor[0]=0;coor[0]< shrd.fout[1].sizes[0];coor[0]++) shrd.fout[1](coor) = exp(bestguess.d[coor[0] + deltastate.sizes[0]+1]);
			if (auto ite = deltastate.getIterator()) do{ *ite = bestguess.d[ite()[0]] - shrd.gpa.hstate(ite());} while(ite++);


			for(coor[0] = 0 ; coor[0] < shrd.nbgenes; coor[0]++){
				shrd.fout[2](coor) = 0.0;
				for(int c =0; c< shrd.nbcelltypes;c++) {
					if (c * shrd.nbgenes + coor[0] > shrd.gpa.gps.getSize()) {printf("outside range %i\n", c * shrd.nbgenes + coor[0] ); break;}
					shrd.fout[2](coor) += shrd.gpa[c * shrd.nbgenes + coor[0]].getMeanAndDeriv_p1(bestguess.d[2], exp(bestguess.d[c+3]),comp_deriv, deltastate, computation_buffer());
				}
			}
			LL += bestguess.k;
			shrd.LL_out[coor[1] ] = bestguess.k;
			shrd.Step_out[coor[1] ] = bestwhich[0];
			shrd.Redo_out[coor[1] ] = bestwhich[1];
		}while(jobite++);
	}
	void operator<<=(const TaskDelinate& other){LL += other.LL; local_deriv += other.local_deriv;}
	//double  operator()(){return LL;}
	double  operator()(const GlobalScope& shrd){KeyElem<double, Tuple<double>> fout; fout.k = LL; fout.d.toMemmove(local_deriv); return LL;}
};
	if (gps.getSize() == 0) return shrd.fout;

	shrd.mincoor.setSize(hstate.sizes[0]);
	shrd.offcoor.setSize(hstate.sizes[0]);
	if (auto daite = hstate.mkIterator()){
		shrd.mincoor[0] = shrd.offcoor[0] = *daite;
		for(uint32_t i=1;i<hstate.sizes[0];i++) {daite++; shrd.mincoor[i] = shrd.offcoor[i] = *daite;}
		while(daite++){
			if (*daite < shrd.mincoor[0]) shrd.mincoor[0] = *daite;
			else if (*daite > shrd.offcoor[0]) shrd.offcoor[0] = *daite;
			for(uint32_t i=1;i<hstate.sizes[0];i++) {daite++;
				if (*daite < shrd.mincoor[i]) shrd.mincoor[i] = *daite;
				else if (*daite > shrd.offcoor[i]) shrd.offcoor[i] = *daite;
			}
		}
	}

	shrd.offcoor -= shrd.mincoor;
	shrd.gridsize =8;
	//shrd.offcoor /= shrd.gridsize - 1;
	//printf("minmax\n");
	shrd.mincoor.show();
	shrd.offcoor.show();
	printf("minmax\n");

	double result = foreach[observ.mkVectorIterator()].run<TaskDelinate>(shrd);
	if (auto ite = shrd.LL_out.mkIterator()) do{
		printf("LL[%i]=%e at(%i,%i)\n", ite(), *ite,  shrd.Step_out[ite()] ,  shrd.Redo_out[ite()] );
	}while(ite++);
return shrd.fout;}

TMatrix<double> GaussianProcessArray::makeGParammatrix()const{ TMatrix<double> fout;
	if (gps.getSize() == 0) return fout;
	fout.setSizes(gps.getSize(), 3 + gps[0].scale.totsize());
	Tuple<uint32_t, 2u> coor;
	for(coor[0] = 0;coor[0] < gps.getSize();coor[0]++){
		coor[1] = 0; fout(coor) = gps[coor[0]].getMean();
		coor[1] = 1; fout(coor) = gps[coor[0]].getNoise();
		coor[1] = 2; fout(coor) = gps[coor[0]].getSignal();
		for(coor[1] = 3; coor[1]<gps[0].scale.totsize()+3;coor[1]++) fout(coor) = gps[coor[0]].scale.data[coor[1]-3];
	}
return fout;}

/*
	void parse(FILE* f){
		int c = fgetc(f);
		Vector< Tuple<char, 4> > scope;
		while(c != EOF){
			switch(c){
				case (int)'<':
					scope.push_back();	fscanf(f,"%.4s", (char*)(scope[scope.size()-1]) ); break;
				break;
			}

			c = fgetc(f);
		}

	}*/

/*
    for(v_in[1]=0;v_in[1]< bone_weights.heap.getSize();v_in[1]++){
        n_in[0] = bone_weights.heap[v_in[1]].first.d[0].d;
        for(v_in[0]=1;v_in[0]<bone_weights.heap[v_in[1]].first.d.getSize();v_in[0]++) n_in[0] += bone_weights.heap[v_in[1]].first.d[v_in[0]].d;
        bone_weights.heap[v_in[1]].first.d[0].d = 255.0f;
        for(v_in[0]--;v_in[0]!=0;v_in[0]--) {
            bone_weights.heap[v_in[1]].first.d[v_in[0]].d *= 255.0f / n_in[0];
            bone_weights.heap[v_in[1]].first.d[0].d -= (float)(uint8_t)bone_weights.heap[v_in[1]].first.d[v_in[0]].d;
        }
    }*/


/*
        // curref *= 2.0f * max_size;

        // ExOp::toOne(curref.ori);

//        ExOp::toOne(curref[0].scale);
//        ExOp::toOne(curref[0].scale);

        curref[0].wrMatrix(dami[0]);
        curref[1].wrMatrix(dami[1]);
        dami[0] *= (1.0f - mix[2]);
        dami[1] *= mix[2];
		dami[0] += dami[1];
		dami[1] = dami[0].mkInverse();


        for(unsigned int c=0;c< v_pos.getSize();c++){
            n_in[0] = v_pos[c][i][0] / (2.0f * max_size);
            n_in[1] = v_pos[c][i][1] / (2.0f * max_size);
            n_in[2] = v_pos[c][i][2] / (2.0f * max_size);
	//n_in[3] =1.0f;
			n_in -= curref[0].pos * (1.0f - mix[2]);
			n_in -= curref[1].pos * mix[2];
			n_in = dami[1] * n_in;



            //v_pos[c][i] -= curref.pos;
    //        v_pos[c][i][0] = n_in[0] * (2.0f * max_size);
  //          v_pos[c][i][1] = n_in[1] * (2.0f * max_size);
   //         v_pos[c][i][2] = n_in[2] * (2.0f * max_size);
    //    }

*/

/*
LockPtr<void, true>::LockPtr(const AliasPtr<void> &tar){initAlias(tar.alias);}
LockPtr<void, true>::LockPtr(const unsigned int &alias){initAlias(alias);}
LockPtr<void, true>::~LockPtr(){if (target != NULL) {
        unsigned int ite = AliasOf.find(target);
        ite = AliasBank.find(AliasOf.deref(ite));
        AliasBank.deref(ite).second -= 0x01000000;}
    }
void LockPtr<void, true>::initAlias(const unsigned int &alias){
   // unsigned int tridmask = ((unsigned int)Controlstate::ThreadID_mask[SDL_ThreadID()]) << 28;
   unsigned int ite = AliasBank.find(alias);
    if (ite == 0xFFFFFFFF) target =NULL;
    else{
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        AliasBank.deref(ite).second += 0x01000000;
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        if ((AliasBank.deref(ite).second & 0xFF000000) == 0x01000000) target = (const void*) AliasBank.deref(ite).first;
        else{AliasBank.deref(ite).second -= 0x01000000;
        target = NULL;}
    }
    }
LockPtr<void,false>::LockPtr(const AliasPtr<void> &tar){initAlias(tar.alias);}
LockPtr<void,false>::LockPtr(const unsigned int &alias){initAlias(alias);}

LockPtr<void,false>::~LockPtr(){if (target != NULL) {
        unsigned int ite = AliasOf.find(target);
        ite = AliasBank.find(AliasOf.deref(ite));
        AliasBank.deref(ite).second -= 0x01000000;
        if ((AliasBank.deref(ite).second & 0x0F000000) == 0) AliasBank.deref(ite).second &= 0x0FFFFFFF;
    }
    }
void LockPtr<void,false>::initAlias(const unsigned int &alias){
   unsigned int tridmask = ((unsigned int)Controlstate::ThreadID_mask[SDL_ThreadID()]) << 28;
   unsigned int ite =  AliasBank.find(alias);// printf("found %i\n", ite);
    if (ite == 0xFFFFFFFF) target =NULL;
    else{
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        if ((AliasBank.deref(ite).second & 0xF0000000) == 0xF0000000) { // read-only lock!
            target = NULL;
        }else if ((AliasBank.deref(ite).second & 0xFF000000) != 0) {
            if  ((AliasBank.deref(ite).second & 0xF0000000) == tridmask) {AliasBank.deref(ite).second += 0x01000000; target = (void*) AliasBank.deref(ite).first;}
            else target = NULL;
        }else{
            AliasBank.deref(ite).second += 0x01000000;
            if ((AliasBank.deref(ite).second & 0xFF000000) == 0x01000000) {target = (void*) AliasBank.deref(ite).first; AliasBank.deref(ite).second |= tridmask;}
            else{ AliasBank.deref(ite).second -= 0x01000000;target = NULL;}
        }
        }
}*/



} // end of namespace

