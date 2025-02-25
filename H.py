import ee
global_ee = ee

def calc_aerodinamic_res(image, wind_image):
    """ 1. AERODINAMIK QARSHILIKNI HISOBLAYMIZ [Shamol tezligi ham kerak bo'ladi] """
    # Ma'lumotlarni olish
    lai = image.select('LAI')
    wind_speed = wind_image.select('u_component_of_wind_10m')  #  TOGRILANDI

    # Z_OM ni hisoblash
    z_om = ee.Image.constant(0.018).multiply(lai)

    # SHAMOL TEZLIGI VA DOIMIY
    wind_height = ee.Image.constant(10)
    k = ee.Image.constant(0.41)  # von Kármán constant

    # Friction velocity ni hisoblash
    u_f_v = k.multiply(wind_speed).divide(wind_height.divide(z_om).log())

    # Z lar qiymatini hisoblash
    z_2 = ee.Image.constant(2)  # 
    z_1 = ee.Image.constant(0.1)

    # Aerodynamic resistance (R_AH) ni hisoblash
    r_ah = ((z_2.divide(z_1)).log()).divide(u_f_v.multiply(k)).rename('R_AH')

    return image.addBands(r_ah)

def calc_different_temp(image, air_t):
        """2. Surface temperature and air temperature lari o'rtasidagi farq [dT ni hisoblash funksiyasi] """
        s_t = image.select('T_S')  # Surface temperature
        a_t = air_t.select('temperature_2m')  # Air temperature 2m dagi
        d_t = s_t.subtract(a_t).rename('dT')  # Haroratlar farqi
        return image.addBands(d_t)  # yangi band qo'shish
    
def calc_air_bosim(image, bosim, harorat):
        """3. Havo zichligini hisoblash [ Meto ma`lumotlar bo'lmasa ishlatiladi ] """
        p = bosim.select('surface_pressure')  # Bosimni tanlash
        t = harorat.select('temperature_2m')  # Haroratni tanlash
        r_d = ee.Image.constant(287.5)  # Gaz konstantasi (havo uchun)

        # Havo zichligini hisoblash
        r = image.expression(
            'p / (R * t)',  # Bu yerda R va t o'zgaruvchilariga to'g'ri ishlov berildi
            {
                'p': p,
                't': t,
                'R': r_d
            }
        ).rename('Bosim')  # Havo zichligini nomlash

        return image.addBands(r)  # R bandini tasvirga qo'shish

def calc_heat_flux(image, *args):
        """4. H ni hisoblash funksiyasi """
        t = image.select('T_S')  # Yuzaning harorati
        r = image.select('R_AH')  # Aerodinamik qarshilik
        cp = ee.Image.constant(1004)  # Havo doimiysi (J/(kg·K))
        dt = image.select('dT')  # Harorat farqi (dT)
        
        if args:
            bosim = args
        else:
            bosim = image.select('Bosim')  # Havo zichligi

        h = image.expression(
            'r * cp * dT / rah',  # Heat flux formula
            {
                'cp': cp,
                'r' : bosim,
                'dT': dt,
                'rah': r
            }
        ).rename('H')  # Natijaviy bandni nomlash

        return image.addBands(h)  # H bandini tasvirga qo'shish


import ee
global_ee = ee

def calc_aerodinamic_res(image, wind_image,*landcover_image):
    """ 1. AERODINAMIK QARSHILIKNI HISOBLAYMIZ [Shamol tezligi ham kerak bo'ladi] /
        aerodynamic resistance to heat transport (rah)                           """
    # rah = ln(z2/z1) / (u_f * k) # rah - aerodynamic resistance to heat transport, u_f - friction velocity, k - von Kármán constant
    # u_f = k * ux / ln(z_x/z_om) # u_f - friction velocity(m/s), ux - wind speed at height z-->z_x, z_om - momentum roughness length(m),
    # Ma'lumotlarni olish
    lai = image.select('LAI')
    wind_speed = wind_image.select('u_component_of_wind_10m')  #  TOGRILANDI
    
    # Land Use - Map dan foydalanib 
    # Z_OM ni hisoblash [ agricultural areas ] uchun
    z_om = ee.Image.constant(0.018).multiply(lai)  # suv, shahar, o'rmon uchun boshqa qiymat bo'ladi 
    
    """ 
    # Har bir sinf uchun maskalarni yaratamiz (Bu Klassifikatsiya qilingan bandlar uchun):
    landcover = landcover_image.select('landcover')  # Bu yerda landcover bandini almashtirish kerak
    water_mask = landcover.eq(1)
    snow_mask = landcover.eq(2)
    desert_mask = landcover.eq(3)
    agriculture_mask = landcover.eq(4)
    bare_soil_mask = landcover.eq(5)
    alfalfa_mask = landcover.eq(6)
    rock_mask = landcover.eq(7)
    
    # Har bir sinf uchun doimiy z_om qiymatlarini aniqlaymiz:
    z_om_water = ee.Image.constant(0.0005)
    z_om_cities = ee.Image.constant(0.2)
    z_om_forest = ee.Image.constant(0.5)  # 
    z_om_grassland = ee.Image.constant(0.02)  # 
    z_om_desert_veg = ee.Image.constant(0.1)
    z_om_snow = ee.Image.constant(0.005)
    
    """
    
    # SHAMOL TEZLIGI VA DOIMIY
    wind_height = ee.Image.constant(10)
    k = ee.Image.constant(0.41)  # von Kármán constant

    # Friction velocity ni hisoblash
    u_f_v = k.multiply(wind_speed).divide(wind_height.divide(z_om).log())

    # Z lar qiymatini hisoblash
    z_2 = ee.Image.constant(2)  # 
    z_1 = ee.Image.constant(0.1)

    # Aerodynamic resistance (R_AH) ni hisoblash
    r_ah = ((z_2.divide(z_1)).log()).divide(u_f_v.multiply(k)).rename('R_AH')

    return image.addBands(r_ah)

def calc_aero_res_with_meto(image, wind_speed, wind_height *args):
    """ 1. AERODINAMIK QARSHILIKNI HISOBLAYMIZ [ Metrologik ma`lumotlar orqali ]"""
    wind_speed = ee.Image.constant(wind_speed)    # speed (m/s)
    lai = image.select('LAI')
    z_om = ee.Image.constant(0.018).multiply(lai)
    wind_height = ee.Image.constant(wind_height)    # wind_height (m)
    k = ee.Image.constant(0.41)  # von Kármán constant
    
    pass # ... to be continued ...------------->>>>>>>>>>>>>>>>>

def calc_different_temp(image, air_t):
        """2. Surface temperature and air temperature lari o'rtasidagi farq [dT ni hisoblash funksiyasi] """
        s_t = image.select('T_S')  # Surface temperature
        a_t = air_t.select('temperature_2m')  # Air temperature 2m dagi
        d_t = s_t.subtract(a_t).rename('dT')  # Haroratlar farqi
        return image.addBands(d_t)  # yangi band qo'shish
    
def calc_air_bosim(image, bosim, harorat):
        """3. Havo zichligini hisoblash [ Meto ma`lumotlar bo'lmasa ishlatiladi ] """
        p = bosim.select('surface_pressure')  # Bosimni tanlash
        t = harorat.select('temperature_2m')  # Haroratni tanlash
        r_d = ee.Image.constant(287.5)  # Gaz konstantasi (havo uchun)

        # Havo zichligini hisoblash
        r = image.expression(
            'p / (R * t)',  # Bu yerda R va t o'zgaruvchilariga to'g'ri ishlov berildi
            {
                'p': p,
                't': t,
                'R': r_d
            }
        ).rename('Bosim')  # Havo zichligini nomlash

        return image.addBands(r)  # R bandini tasvirga qo'shish

def calc_heat_flux(image, *args):
        """4. H ni hisoblash funksiyasi """
        t = image.select('T_S')  # Yuzaning harorati
        r = image.select('R_AH')  # Aerodinamik qarshilik
        cp = ee.Image.constant(1004)  # Havo doimiysi (J/(kg·K))
        dt = image.select('dT')  # Harorat farqi (dT)
        
        if args:
            bosim = args
        else:
            bosim = image.select('Bosim')  # Havo zichligi

        h = image.expression(
            'r * cp * dT / rah',  # Heat flux formula
            {
                'cp': cp,
                'r' : bosim,
                'dT': dt,
                'rah': r
            }
        ).rename('H')  # Natijaviy bandni nomlash

        return image.addBands(h)  # H bandini tasvirga qo'shish

"""
SEBAL computes dT for each pixel by assuming a linear relationship between dT and Ts:
                              dT = b + aTs
# b and a are the correlation coefficients
"""
