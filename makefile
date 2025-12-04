# 定义编译器
CXX = mpicxx

# 定义编译选项, 编译步骤通常使用CXXFLAGS，而链接步骤使用LDFLAGS
CXXFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 -march=native -O3
INCLUDES = -I ../ccp4_lib/ 

# 定义链接选项 #/usr/local/lib/libfftw3.dylib, -lfftw3已经包括了-lfftw和-lrfftw
LDFLAGS = -L ../ccp4_lib/lib -L /usr/local/fftw2/lib
#LIBS =  -lclipper-ccp4 -lclipper-core -lclipper-contrib -lccp4c -lccp4f -lfftw3
LIBS =  -lclipper-ccp4 -lclipper-core -lclipper-contrib -lccp4c -lccp4f -lfftw -lrfftw

# 目标文件和源文件路径
SRC_DIR = src
OBJ_DIR = obj
OUTPUT_DIR = output
INCLUDE_DIR = include

# 定义源文件和目标文件
SRC = $(SRC_DIR)/main.cpp $(SRC_DIR)/functions.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/globalvars.cpp $(SRC_DIR)/dataimporter.cpp
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))
TARGET = program

# 规则：生成最终的可执行文件，@代表目标文件的名称，也就是规则中的(TARGET)部分。^代表所有的依赖项，也就是(OBJ)列表中的所有文件。
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

# 编译每个 .cpp 文件为 .o 文件
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# 创建 obj, output 目录
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR) $(OUTPUT_DIR)

# 清理生成的文件
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(TARGET)
