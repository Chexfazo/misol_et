try:
    import geemap
    import ee
    
except ModuleNotFoundError:
    st.error("geemap kutubxonasini o‘rnatib bo‘lmadi. Iltimos, requirements.txt faylni tekshiring!")

def calc_spek_rad(image):
    """Landsat 8 tasviri uchun spektral radiance hisoblash funksiyasi"""

    for x in range(1, 12):  # B1 dan B11 gacha  0 dan boshlanmaydi !
        qcal = image.select(f'B{x}')
        m_l = ee.Number(image.get(f'RADIANCE_MULT_BAND_{x}'))
        add = ee.Number(image.get(f'RADIANCE_ADD_BAND_{x}'))

        # Spektral radiance hisoblash
        radiance = qcal.multiply(m_l).add(add)

        # Yangi band qo‘shish
        image = image.addBands(radiance.rename(f'SR_{x}'))

    return image

def reflectivity(image):
        """Landsat 8 tasviri uchun Reflectivityni hisoblash funksiyasi"""
        # pi qiymati
        pi = ee.Number(3.141592)

        # dr ni hisoblash
        doy = ee.Number(ee.Date(image.get('system:time_start')).getRelative('day', 'year')).add(1)
        dr = ee.Number.expression(
            '1 + 0.033 * cos(DOY * 2 * PI / 365)',
            {
                'DOY': doy,
                'PI': pi
            }
        )

        # cos() ni hisoblash
        sun_elev = ee.Number(image.get('SUN_ELEVATION'))
        alfa = ee.Number(90).subtract(sun_elev)
        cos_alfa = alfa.multiply(pi).divide(180).cos()

        # esun qiymatlari
        esun_org = [1895.33, 2004.57, 1820.75, 1549.49, 951.76, 247.55, 85.46, 1723.88, 366.97, 1, 1]  # list 0 dan boshlanadi , dict kabi yozsa ham bo`ladi{1:1890, 2:2004, }

        # LOOP yozish
        for n in range(1, 12):  # 1 dan 11 gacha bandlar
            esun = ee.Number(esun_org[n-1])  # esun ni olish natija 0 sababi list 0 dan boshlandi
            band = image.select(f'SR_{n}')

            # Reflectivity ni hisoblash
            ref = image.expression(
                '(pi * band) / (esun * dr * cos_alfa)',
                {
                    'pi': pi,
                    'band': band,
                    'esun': esun,
                    'dr': dr,
                    'cos_alfa': cos_alfa
                }
            )

            #   bandni Nomi bilan qo‘shish
            image = image.addBands(ref.rename(f'REF_{n}'))

        return image
    # yakuniy natija REF_1, REF_2 ... kabi boladi
    
def albedo_toa(image):
    """Landsat 8 tasviri uchun Albedo at the Top Atmosphere (TOA) ni hisoblash funksiyasi"""

    # Landsat 8 uchun ESUN qiymatlari (B1-B11)
    esun_org = [1895.33, 2004.57, 1820.75, 1549.49, 951.76, 247.55, 85.46, 1723.88, 366.97, 1, 1]

    # jami ESUN ni hisoblash
    sum_esun = sum(esun_org[:10])  # Faqat B1-B10 ni hisobga olamiz, B11 = 1

    # Albedo TOA hisoblash uchun boshlang'ich qiymatlar
    weighted_sum = ee.Image(0)

    # Barcha bandlar bo‘yicha sikl
    for n in range(1, 11):  # B1-B10 uchun ishlaydi, B11 termal band (kerak emas)
        esun = ee.Number(esun_org[n - 1])  # `n - 1` qilish kerak, chunki Python indekslash 0 dan boshlanadi
        w = esun.divide(sum_esun)  # W_i = ESUN_i / ∑ESUN

        # REF bandni olish
        ref_band = image.select(f'REF_{n}')

        # TOA Albedo uchun hisoblash
        weighted_sum = weighted_sum.add(ref_band.multiply(w))

        # Yakuniy natija
        albedo_toa = weighted_sum.rename('ALBEDO_TOA')

    return image.addBands(albedo_toa)

def calc_albedo(image, srtm):
        """Landsat 8 tasviri uchun Albedo ni hisoblash funksiyasi"""
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        tsw = tsw.pow(2)
        path_radian = ee.Image.constant(0.03)  # Bastiaansen taklif qilgan qiymat
        toa = image.select('ALBEDO_TOA')
        albedo = (toa.subtract(path_radian)).divide(tsw)

        return image.addBands(albedo.rename('ALBEDO'))
    
def calc_inshort_radiation(image, srtm):
        """Landsat 8 tasviri uchun Incoming Shortwave Radiationni hisoblash funksiyasi"""
        pi = ee.Number(3.141592)

        # tsw hisoblash
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        tsw = tsw.pow(2)

        # dr - Earth-Sun distance correction factor (IMAGE ga o‘giramiz)
        doy = ee.Number(ee.Date(image.get('system:time_start')).getRelative('day', 'year')).add(1)
        dr = ee.Image.constant(ee.Number(1).add(ee.Number(0.033).multiply(doy.multiply(ee.Number(2)).multiply(pi).divide(365).cos())))

        # cos() hisoblash
        sun_elev = ee.Number(image.get('SUN_ELEVATION'))
        alfa = ee.Number(90).subtract(sun_elev)
        cos_alfa = ee.Image.constant(alfa.multiply(pi).divide(180).cos())

        # G (quyosh doimiysi) IMAGE sifatida
        g = ee.Image.constant(1367)

        # Yakuniy natija IMAGE bo‘lishi kerak
        inshort_radiation = g.multiply(dr).multiply(cos_alfa).multiply(tsw)

        return image.addBands(inshort_radiation.rename('IN_SHORT_RAD'))
    
def calc_surface_ems(image):
        """Landsat 8 tasviri uchun Surface Emissivityni hisoblash funksiyasi"""

        nir = image.select('REF_5')
        red = image.select('REF_4')

        # NDVI hisoblash
        ndvi = image.normalizedDifference(['REF_5', 'REF_4']).rename('NDVI')

        # SAVI hisoblash
        l = ee.Image.constant(0.5)
        savi = image.expression(
            '1.5 * (nir - red) / (nir + red + l)',
            {
                'nir': nir,
                'red': red,
                'l': l
            }
        ).rename('SAVI')

        # LAI hisoblash
        lai = ee.Image.constant(0.69).subtract(savi).divide(0.59).log().divide(ee.Image.constant(0.91)).rename('LAI')

        # Surface emissivity E_NB va E_O
        e_nb = ee.Image(0.98).where(ndvi.gt(0).And(lai.lt(3)), ee.Image(0.97).add(ee.Image.constant(0.0033).multiply(lai))).rename('E_NB')
        e_o = ee.Image(0.98).where(ndvi.gt(0).And(lai.lt(3)), ee.Image(0.95).add(ee.Image.constant(0.01).multiply(lai))).rename('E_O')

        # Natijalarni tasvirga qo'shish
        return image.addBands([e_nb, e_o, lai, ndvi, savi])
    
def calc_surface_temp(image):
        """Landsat 8 tasviri uchun Surface Temperature hisoblash funksiyasi"""

        # K1 va K2 konstantalarni olish
        k_1 = ee.Number(image.get('K1_CONSTANT_BAND_10'))
        k_2 = ee.Number(image.get('K2_CONSTANT_BAND_10'))

        # Temperature bandini olish
        radiance = image.select('REF_10')

        # Surface Temperature (T_S) hisoblash
        t_s = image.expression(
            'K2 / log(ENB * K1 / Rc + 1)',
            {
                'K2': k_2,
                'ENB': image.select('E_NB'),
                'K1': k_1,
                'Rc': radiance
            }
        ).rename('T_S')

        # Natijani image ga qo'shish
        return image.addBands(t_s)
    
def cal_out_rad(image):
        """Landsat 8 tasviri uchun OUTGOING LONGWAVE RADIATION (OLR) hisoblash funksiyasi"""

        # Kerakli bandlarni tanlash
        t_s = image.select('T_S')  # Surface Temperature (Kelvin)
        e_o = image.select('E_O')  # Surface emissivity

        # Stefan-Boltzmann doimiysi (W/m^2K^4)
        stef_bolt = ee.Image.constant(5.67e-8)

        # Outgoing Longwave Radiation (OLR) hisoblash
        out_rad = stef_bolt.multiply(e_o).multiply(t_s.pow(4)).rename('OUT_RAD')

        # Natijani image ga qo‘shish
        return image.addBands(out_rad)

def calc_in_long(image, air_temp, srtm):
        """Landsat 8 tasviri uchun Incoming LONGWAVE RADIATION (OLR) hisoblash funksiyasi"""

        # tsw (atmosferik transmittans) ni hisoblash
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        tsw = tsw.pow(2)

        # e_a (atmosferik emissivlik) ni hisoblash
        e_a = ee.Image.constant(0.85).multiply(tsw.pow(0.09).log().multiply(-1))

        # Stefan-Boltzmann doimiysi
        step_boltz = ee.Image.constant(5.67e-8)

        # Havo harorati (air_temp ni tanlash)
        t_air = air_temp.select('temperature_2m')

        # Incoming Longwave Radiation (OLR) hisoblash
        in_long = step_boltz.multiply(e_a).multiply(t_air.pow(4)).rename('IN_LONG')

        # Natijani image ga qo'shish
        return image.addBands(in_long)
    
def net_radiation(image):
        """Landsat 8 tasviri uchun NET RADIATION NI(Rn) hisoblash funksiyasi"""
        albedo = image.select('ALBEDO')
        inshort_rad = image.select('IN_SHORT_RAD')
        out_rad = image.select('OUT_RAD')
        in_long = image.select('IN_LONG')
        e_o = image.select('E_O')

        # Rn ni hisoblash
        net_rad = image.expression(
            '(1-a)*RS - RL_T + RL_P * E_O',
            {
                'a': albedo,
                'RS': inshort_rad,
                'RL_T': out_rad,
                'RL_P': in_long,
                'E_O': e_o
            }
        ).rename('Rn')

        return image.addBands(net_rad)
    
def soil_heat_flux(image):
        """ G ni SOIL HEAT flux ni hisoblaymiz """
        r_n = image.select('Rn')
        t_s = image.select('T_S')
        a = image.select('ALBEDO')
        ndvi = image.select('NDVI')
        ndvi = ndvi.pow(4)

        g = image.expression(
            '(Rn * Ts) / a * (0.0038 * a + 0.0074 * a * a) * (1 - 0.98 * ndvi)',
            {
                'Rn': r_n,
                'Ts': t_s,
                'a': a,
                'ndvi': ndvi
            }
        ).rename('G')

        # band qo'shish
        return image.addBands(g)
    
def calc_aerodinamic_res(image, wind_image):
    """ AERODINAMIK QARSHILIKNI HISOBLAYMIZ """
    # Ma'lumotlarni olish
    lai = image.select('LAI')
    wind_speed = wind_image.select('u_component_of_wind_10m')  # ✅ TO‘G‘RILANDI

    # Z_OM ni hisoblash
    z_om = ee.Image.constant(0.018).multiply(lai)

    # SHAMOL TEZLIGI VA DOIMIY
    wind_height = ee.Image.constant(10)
    k = ee.Image.constant(0.41)  # von Kármán constant

    # Friction velocity ni hisoblash
    u_f_v = k.multiply(wind_speed).divide(wind_height.divide(z_om).log())

    # Z lar qiymatini hisoblash
    z_2 = ee.Image.constant(2)  # ✅ CONSTANTGA O‘GIRILDI
    z_1 = ee.Image.constant(1)

    # Aerodynamic resistance (R_AH) ni hisoblash
    r_ah = (z_2.log().subtract(z_1.log())).divide(u_f_v.multiply(k)).rename('R_AH')

    return image.addBands(r_ah)


def different_temp(image, air_t):
        """ Surface temperature and air temperature lari o'rtasidagi farq"""
        s_t = image.select('T_S')  # Surface temperature
        a_t = air_t.select('temperature_2m')  # Air temperature 2m dagi
        d_t = s_t.subtract(a_t).rename('dT')  # Haroratlar farqi
        return image.addBands(d_t)  # yangi band qo'shish

def calc_havo_bosim(image, bosim, harorat):
        """ Havo zichligini hisoblash """
        p = bosim.select('mean_sea_level_pressure')  # Bosimni tanlash
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

def calc_heat_flux(image):
        """ H ni hisoblash funksiyasi """
        t = image.select('T_S')  # Yuzaning harorati
        r = image.select('R_AH')  # Aerodinamik qarshilik
        cp = ee.Image.constant(1004)  # Havo doimiysi (J/(kg·K))
        dt = image.select('dT')  # Harorat farqi (dT)

        h = image.expression(
            'cp * dT / rah',  # Heat flux formula
            {
                'cp': cp,
                'dT': dt,
                'rah': r
            }
        ).rename('H')  # Natijaviy bandni nomlash

        return image.addBands(h)  # H bandini tasvirga qo'shish

def evapotranspiration(image):
        """ Evapotranspiratsiyani hisoblash funksiyasi """
        bands_needed = ['Rn', 'G', 'H']

        # Tasvirda kerakli bandlar borligini tekshiramiz
        missing_bands = [band for band in bands_needed if not image.bandNames().contains(band)]

        if missing_bands:
            raise ValueError(f"Quyidagi bandlar yetishmayapti: {missing_bands}")

        r_n = image.select('Rn')
        g = image.select('G')
        h = image.select('H')

        # ET = Rn - G - H
        et = r_n.subtract(g).subtract(h).rename('ET')

        return image.addBands(et)



Map = geemap.Map()
Map

landsat = ee.Image('LC08_156033_20190317').filterBounds(roi)
roi = Map.draw_features
roi = ee.FeatureCollection(roi)

air_temp = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY").filterBounds(roi).filterDate('2019-01-01', '2019-12-31').first()
  
wind_image = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY") \
    .filterBounds(roi) \
    .filterDate('2019-01-01', '2019-12-31') \
    .first()  # Birinchi tasvirni tanlaymiz

srtm = ee.Image("USGS/SRTMGL1_003").clip(roi)
  
bosim = ee.ImageCollection("ECMWF/ERA5/DAILY") \
    .filterBounds(roi) \
    .filterDate('2019-01-01', '2019-12-31').first()


def result_final(image,srtm,air_temp,wind_image,bosim):
    """Yakuniy natijani hisoblaymiz Azizlar !"""
    spek_rad = calc_spek_rad(image) 
    reflictivity = reflectivity(spek_rad)
    top_albedo = albedo_toa(reflictivity)
    albedo = calc_albedo(top_albedo)
    in_short = calc_inshort_radiation(albedo)
    emissivity = calc_surface_ems(in_short)
    surface_temp = calc_surface_temp(emissivity)
    out_rad = cal_out_rad(surface_temp)
    in_long = calc_in_long(out_rad,air_temp, srtm)
    net_radian = net_radiation(in_long)
    soil_heat = soil_heat_flux(net_radian)
    aero_res = calc_aerodinamic_res(soil_heat,wind_image)
    different_t = different_temp(aero_res, air_temp)
    air_pressure = calc_havo_bosim(different_t, bosim, wind_image)
    heat_flux = calc_heat_flux(air_pressure)
    e_t = evapotranspiration(heat_flux)

result_final(landsat,srtm,air_temp,wind_image,bosim)
Map.addLayer(landsat, {}, 'LANDSAT')
Map

