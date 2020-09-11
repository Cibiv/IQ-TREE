add_library(terrace
terracenode.cpp
terracenode.hpp
terracetree.cpp
terracetree.hpp
presenceabsencematrix.cpp
presenceabsencematrix.hpp
)

target_link_libraries(terrace tree)
