wasm:
	emcc -O3 ../geokit.c -o geokit.js \
		-sALLOW_MEMORY_GROWTH=1 \
		-sWASM_ASYNC_COMPILATION=1 \
		-sMODULARIZE \
		-sEXPORTED_RUNTIME_METHODS=ccall \
		-sEXPORTED_FUNCTIONS="[\
			'_geokit_inverse',\
			'_geokit_direct',\
			'_geokit_track',\
			'_geokit_intercept',\
			'_geokit_intersect',\
			'_geokit_nearest',\
			'_geokit_equibox',\
			'_geokit_equipath',\
			'_malloc',\
			'_free'\
		]"
