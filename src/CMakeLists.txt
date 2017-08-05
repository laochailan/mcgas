set(SRCs
	mcgas.c
	random.c
	render.c
)

find_package(OpenMP REQUIRED)
find_package(WebP REQUIRED)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3 -ffast-math ${OpenMP_C_FLAGS}")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -O3 -ffast-math")

include_directories(${WEBP_INCLUDE_DIRS})
add_executable(mcgas ${SRCs})


target_link_libraries(mcgas -lm ${WEBP_LIBRARIES})
