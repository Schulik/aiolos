#ifndef _AIOLOS_ENUM_H_
#define _AIOLOS_ENUM_H_

#include <fstream>

#define IOS_INPUT(ENUM_TYPE) \
    inline std::istream& operator>>(std::istream& is, ENUM_TYPE& obj) { \
        std::underlying_type<ENUM_TYPE>::type tmp ; \
        is >> tmp ; \
        obj = static_cast<ENUM_TYPE>(tmp) ; \
        return is ; \
    }

#define IOS_OUTPUT(ENUM_TYPE) \
    inline std::ostream& operator<<(std::ostream& os, ENUM_TYPE obj) { \
        os << static_cast<std::underlying_type<ENUM_TYPE>::type>(obj) ;  \
        return os ; \
    }

enum class Geometry {
    cartesian = 0, cylindrical = 1, spherical = 2
} ;
IOS_INPUT(Geometry);
IOS_OUTPUT(Geometry);


enum class BoundaryType {
    user = 0, open = 1, reflecting = 2, fixed = 3, periodic = 4,
} ;
IOS_INPUT(BoundaryType);
IOS_OUTPUT(BoundaryType);

enum class IntegrationType {
    first_order = 1, second_order = 2, WENO=3
} ;
IOS_INPUT(IntegrationType);
IOS_OUTPUT(IntegrationType);


enum class EOS_pressure_type {
    adiabatic = 0, polytropic = 1, tabulated = 2, supernova = 3, user=4
};
IOS_INPUT(EOS_pressure_type);
IOS_OUTPUT(EOS_pressure_type);

enum class EOS_internal_energy_type {
    thermal = 0, constant = 1, tabulated = 2, supernova = 3, user=4
};
IOS_INPUT(EOS_internal_energy_type);
IOS_OUTPUT(EOS_internal_energy_type);



#endif//_AIOLOS_ENUM_H