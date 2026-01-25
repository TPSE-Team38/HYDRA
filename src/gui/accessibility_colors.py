GreenBox="#c6f6d5"
Normal_GreenBox="#c6f6d5"
Deuteranomaly_GreenBox="#cfeef2"
Deuteranopia_GreenBox="#bfe6ee"
Protanomaly_GreenBox="#d4eef2"
Protanopia_GreenBox="#b6dde6"
Tritanomaly_GreenBox="#e6f2cc"
Tritanopia_GreenBox="#d9efb3"
Cone_Monochromacy_GreenBox="#dfdfdf"
Achromatopsia_GreenBox="#d9d9d9"

GreenBox_Border="#2e7d32"
Normal_GreenBox_Border="#2e7d32"
Deuteranomaly_GreenBox_Border="#1f6f7a"
Deuteranopia_GreenBox_Border="#155e6a"
Protanomaly_GreenBox_Border="#2c6f7a"
Protanopia_GreenBox_Border="#1f5f6a"
Tritanomaly_GreenBox_Border="#6a7f2e"
Tritanopia_GreenBox_Border="#5f7324"
Cone_Monochromacy_GreenBox_Border="#5a5a5a"
Achromatopsia_GreenBox_Border="#1a1a1a"


YellowBox="#fdf0ac"
Normal_YellowBox="#fdf0ac"
Deuteranomaly_YellowBox="#fbe9c7"
Deuteranopia_YellowBox="#f5e6b8"
Protanomaly_YellowBox="#f9ebcd"
Protanopia_YellowBox="#efe0b0"
Tritanomaly_YellowBox="#fdeae3"
Tritanopia_YellowBox="#f0e4cc"
Cone_Monochromacy_YellowBox="#eeeeee"
Achromatopsia_YellowBox="#e0e0e0"

YellowBox_Border="#ba8e23"
Normal_YellowBox_Border="#ba8e23"
Deuteranomaly_YellowBox_Border="#a69324"
Deuteranopia_YellowBox_Border="#7a6f2a"
Protanomaly_YellowBox_Border="#a19124"
Protanopia_YellowBox_Border="#6f6524"
Tritanomaly_YellowBox_Border="#ba884a"
Tritanopia_YellowBox_Border="#8a6f4a"
Cone_Monochromacy_YellowBox_Border="#8c8c8c"
Achromatopsia_YellowBox_Border="#333333"


RedBox="#fed7d7"
Normal_RedBox="#fed7d7"
Deuteranomaly_RedBox="#e6d8f0"
Deuteranopia_RedBox="#d8ccee"
Protanomaly_RedBox="#e8d6c8"
Protanopia_RedBox="#dcc4b0"
Tritanomaly_RedBox="#fecfd1"
Tritanopia_RedBox="#fed0d0"
Cone_Monochromacy_RedBox="#e9e9e9"
Achromatopsia_RedBox="#cccccc"

RedBox_Border="#c62828"
Normal_RedBox_Border="#c62828"
Deuteranomaly_RedBox_Border="#6a3d7c"
Deuteranopia_RedBox_Border="#55306a"
Protanomaly_RedBox_Border="#7a4a2e"
Protanopia_RedBox_Border="#6a3a1e"
Tritanomaly_RedBox_Border="#c6284a"
Tritanopia_RedBox_Border="#c62860"
Cone_Monochromacy_RedBox_Border="#4b4b4b"
Achromatopsia_RedBox_Border="#000000"

current_mode="Normal"

def change_accessibility_color(mode:str):
    global GreenBox,GreenBox_Border,YellowBox,YellowBox_Border, RedBox, RedBox_Border,current_mode
    GreenBox=globals()[mode+"_GreenBox"]
    GreenBox_Border=globals()[mode+"_GreenBox_Border"]
    YellowBox=globals()[mode+"_YellowBox"]
    YellowBox_Border=globals()[mode+"_YellowBox_Border"]
    RedBox=globals()[mode+"_RedBox"]
    RedBox_Border=globals()[mode+"_RedBox_Border"]
    current_mode=mode
