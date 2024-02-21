#include <vector>

#include "TColor.h"

namespace JJColor
{
    static const TColor greenWut(TColor::GetFreeColorIndex(),0.,104./255,114./255);
    static const TColor navyWut(TColor::GetFreeColorIndex(),0.,74./255,108./255);
    static const TColor beigeWut(TColor::GetFreeColorIndex(),200./255,193./255,182./255);
    static const std::vector<Color_t> fWutMainColors{
        static_cast<Color_t>(TColor::GetColor(greenWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(navyWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(beigeWut.AsHexString()))
    };

    static const TColor mintWut(TColor::GetFreeColorIndex(),192./255,209./255,200./255);
    static const TColor blueWut(TColor::GetFreeColorIndex(),52./255,102./255,175./255);
    static const TColor goldWut(TColor::GetFreeColorIndex(),231./255,162./255,23./255);
    static const TColor fireRedWut(TColor::GetFreeColorIndex(),230./255,51./255,18./255);
    static const TColor wineRedWut(TColor::GetFreeColorIndex(),154./255,15./255,6./255);
    static const TColor violetWut(TColor::GetFreeColorIndex(),94./255,49./255,86./255);
    static const std::vector<Color_t> fWutSecondaryColors{
        static_cast<Color_t>(TColor::GetColor(mintWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(blueWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(goldWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(fireRedWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(wineRedWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(violetWut.AsHexString()))
    };

    static const std::vector<Color_t> fWutAllColors{
        static_cast<Color_t>(TColor::GetColor(greenWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(navyWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(beigeWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(mintWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(blueWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(goldWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(fireRedWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(wineRedWut.AsHexString())),
        static_cast<Color_t>(TColor::GetColor(violetWut.AsHexString()))
    };

    void CreatePrimaryWutGradient()
    {
        int MyPalette[100];
        double Red[]    = {0.,200./255,0.};
        double Green[]  = {74./255,193./255,104./255};
        double Blue[]   = {108./255,182./255,114./255};
        double Length[] = {0.,0.5,1.};
        // navyWut -> beigeWut -> greenWut
        int fi = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
        for (int i = 0; i < 100; ++i) MyPalette[i] = fi + i;
    }

    void CreateSecondaryWutGradient()
    {
        int MyPalette[100];
        double Red[]    = {192./255,122./255,52./255,142./255,231./255,231./255,230./255,192./255,154./255,124./255,94./255};
        double Green[]  = {209./255,156./255,102./255,132/255,162./255,107./255,51./255,33./255,15./255,32./255,49./255};
        double Blue[]   = {200./255,188./255,175./255,99./255,23./255,21./255,18./255,6./255,6./255,46./255,86./255};
        double Length[] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.};
        // mintWut -> intermidiate -> blueWut -> intermidiate -> goldWut -> intermidiate -> fireRedWut -> intermidiate -> wineRedWut -> intermidiate -> violetWut
        int fi = TColor::CreateGradientColorTable(11,Length,Red,Green,Blue,100);
        for (int i = 0; i < 100; ++i) MyPalette[i] = fi + i;
    }

} // namespace JJColor
