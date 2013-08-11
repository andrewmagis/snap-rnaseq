//
// Lookup tables that are too annoying to initialize with array expressions.
// These includes things like complements of bases, numerical values, etc.
//
// To avoid having an init() function that everyone must call, we stow these
// in a special class and have one static instance of that class and variables
// pointing to it to ensure that the initializer (class constructor) is run.
//

#pragma once


class Tables
{
    char complement[256];
    char isN[256];
    int baseValue[256];
    int baseValueNoN[256];  // Same as above but N maps to 0 instead of 4
    char valueBase[5];

    unsigned isLowerCase[256];
    char toUpperCase[256];

    char packedBaseValue[256];
    char packedQualityMask[256];
    char packedValueBase[256];
    char packedValueBaseRC[256];

public:
    Tables();

    const char *getComplement() const { return complement; }
    const char *getIsN() const { return isN; }
    const int  *getBaseValue() const { return baseValue; }
    const int  *getBaseValueNoN() const { return baseValueNoN; }
    const char *getValueBase() const { return valueBase; }

    const char* getPackedBaseValue() const { return packedBaseValue; }
    const char* getPackedQualityMask() const { return packedQualityMask; }
    const char* getPackedValueBase() const { return packedValueBase; }
    const char* getPackedValueBaseRC() const { return packedValueBaseRC; }
    const unsigned *getIsLowerCase() const {return isLowerCase; }
    const char *getToUpperCase() const {return toUpperCase; }
};

extern const char *COMPLEMENT;
extern const char *IS_N;
extern const int  *BASE_VALUE;
extern const char *VALUE_BASE;
extern const char *PACKED_BASE_VALUE;
extern const char *PACKED_QUALITY_MASK;
extern const char *PACKED_VALUE_BASE;
extern const char *PACKED_VALUE_BASE_RC;
extern const int  *BASE_VALUE_NO_N;
extern const unsigned *IS_LOWER_CASE;
extern const char *TO_UPPER_CASE;
