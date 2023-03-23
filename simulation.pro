CONFIG += qt
#CONFIG += debug
CONFIG += c++11

INCDIR = $$PWD/install/include # Third party header files
LIBDIR = $$PWD/install/lib  # Third party libraries
BINDIR = .

DESTDIR = $$BINDIR
TARGET = model

INCLUDEPATH += $${INCDIR}

QMAKE_CXXFLAGS += -I$${INCDIR}, -DNO_FREETYPE

#Only neccesary for profiling
#QMAKE_CXXFLAGS+=-pg
#QMAKE_CXXFLAGS += -g 
#QMAKE_LFLAGS+=-pg

LIBS += -L$${LIBDIR} -lsbml -lglpk -lpngwriter -lfreetype -lpng

HEADERS = BactParam.h \
          project.h \
          Grid.h \
          Cell.h \
          parameter.h \
          Gut_Output.h \
          Graphics.h

SOURCES = AGORAparam.cpp \
          BP_functions.cpp \
          simulation.cpp \
          Grid.cpp \
          parameter.cpp \
          mixedfunctions.cpp \
          Cell.cpp \
          Gut_Output.cpp \
          Graphics.cpp
          

# finished
