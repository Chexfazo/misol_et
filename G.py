import ee
global_ee = ee

def calc_soil_heat_flux(image):
        """ G ni SOIL HEAT flux ni hisoblaymiz """
        r_n = image.select('Rn')
        t_s = image.select('T_S')
        a = image.select('ALBEDO')
        ndvi = image.select('NDVI')
        ndvi = ndvi.pow(4)

        g = image.expression(
            'sRn * Ts / a * (0.0038 * a + 0.0074 * a * a) * (1 - 0.98 * ndvi)',
            {
                'Rn': r_n,
                'Ts': t_s,
                'a': a,
                'ndvi': ndvi
            }
        ).rename('G')

        # band qo'shish
        return image.addBands(g)