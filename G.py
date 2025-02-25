import ee
global_ee = ee

# To'g`ridan to'g`ri G ni hisoblaydi hech qandya shartlarni e`tibor qilmasdan`
def calc_soil_heat_flux(image):
        """ G ni SOIL HEAT flux ni hisoblaymiz """
        r_n = image.select('Rn')
        t_s = image.select('T_S')
        a = image.select('ALBEDO')
        ndvi = image.select('NDVI')
        ndvi = ndvi.pow(4)

        g = image.expression(
            'Rn * Ts / a * (0.0038 * a + 0.0074 * a * a) * (1 - 0.98 * ndvi)',
            {
                'Rn': r_n,
                'Ts': t_s,
                'a': a,
                'ndvi': ndvi
            }
        ).rename('G')

        # band qo'shish
        return image.addBands(g)
    
# G ni hisoblash barcha shartlarni qo'llab quvatlaydi !  
def calc_GRn_ratio(image):
    """G/Rn, ya'ni tuproq issiqlik oqimi nisbati (G) ni hisoblaymiz."""
    r_n = image.select('Rn')
    t_s = image.select('T_S')  # Surface temperature in Kelvin
    a = image.select('ALBEDO')
    ndvi = image.select('NDVI')
    ndvi_pow = ndvi.pow(4)
    ts_s = t_s.subtract(273.15)  # Celsius ga o'tkazish

    # Asosiy ifoda: (bu yerda o'z formulaingizni aniqlik kiritish zarur)
    g_rn = image.expression(
        'Ts_s / a * (0.0038 + 0.0074 * a) * (1 - 0.98 * ndvi_pow)',
        {
            'Ts_s': ts_s,
            'a': a,
            'ndvi_pow': ndvi_pow
        }
    ).rename('G_Rn')
    
    # Agar NDVI < 0 bo'lsa, yoki T_S (Celsius) < 4 va albedo > 0.45 bo'lsa, G/Rn = 0.5
    g_rn = g_rn.where(ndvi.lt(0), 0.5)
    g_rn = g_rn.where(ts_s.lt(4).And(a.gt(0.45)), 0.5)
    
    # Hisoblangan G (so'nggi qiymat: G = G/Rn * Rn)
    g = g_rn.multiply(r_n).rename('G')
    
    # Natijaviy tasvirga qo'shish
    return image.addBands(g_rn).addBands(g)

def calc_G_Rn_value(image):
    """Turli Yer Qoplamalari Uchun G/Rn Qiymatlari"""
    # Faraz qilaylik, 'landcover' bandida yer qoplamasi sinfi raqamlari mavjud:
    # 1 - Deep, clear water; 2 - Snow; 3 - Desert; 4 - Agriculture; 5 - Bare soil; 6 - Full cover alfalfa; 7 - Rock

    # Har bir sinf uchun doimiy G/Rn qiymatlarini aniqlaymiz:
    g_rn_water = ee.Image.constant(0.5)
    g_rn_snow = ee.Image.constant(0.5)
    g_rn_desert = ee.Image.constant(0.3)  # misol uchun o'rtacha qiymat 0.3
    g_rn_agriculture = ee.Image.constant(0.1)  # misol uchun 0.1
    g_rn_bare_soil = ee.Image.constant(0.3)
    g_rn_alfalfa = ee.Image.constant(0.04)
    g_rn_rock = ee.Image.constant(0.3)

    # Har bir sinf uchun maskalarni yaratamiz:
    landcover = image.select('landcover')  # Bu yerda landcover bandini almashtiring
    water_mask = landcover.eq(1)
    snow_mask = landcover.eq(2)
    desert_mask = landcover.eq(3)
    agriculture_mask = landcover.eq(4)
    bare_soil_mask = landcover.eq(5)
    alfalfa_mask = landcover.eq(6)
    rock_mask = landcover.eq(7)

    # Har bir sinf uchun G/Rn bandini yaratamiz:
    g_rn_image = (g_rn_water.updateMask(water_mask)
                .blend(g_rn_snow.updateMask(snow_mask))
                .blend(g_rn_desert.updateMask(desert_mask))
                .blend(g_rn_agriculture.updateMask(agriculture_mask))
                .blend(g_rn_bare_soil.updateMask(bare_soil_mask))
                .blend(g_rn_alfalfa.updateMask(alfalfa_mask))
                .blend(g_rn_rock.updateMask(rock_mask))
                ).rename('G_Rn_constant')

    # Keyin, asosiy image ga bu bandni qo'shamiz:
    image = image.addBands(g_rn_image)
    
    
def calc_g_for_water_sur(image):
    


# ... to be continued ...------------->>>>>>>>>>>>>>>>>
    """ G ni WATER SURFACE uchun hisoblaymiz """
    # From july to December the instanteous G during midday 
    r_n = image.select('Rn')
    a = ee.Image.constant(90)
    b= ee.Image.constant(100)
    
    g_jul_dec = r_n.subtract(a)
    g_24_hour_jul_dec = r_n.multiply(b)
    
    # Janury to June the instanteous G during midday
    c = ee.Image.constant(40)
    d = ee.Image.constant(50)
    g_jan_jun = r_n.subtract(c)
    g_24_hour_jan_jun = r_n.multiply(d)
    
    return image.addBands(g_jul_dec).addBands(g_24_hour_jul_dec).addBands(g_jan_jun).addBands(g_24_hour_jan_jun)
