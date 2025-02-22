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
