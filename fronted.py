import streamlit as st 
import geemap
import ee
import backgroundkod as bk 
from datetime import datetime

# Streamlit sahifa sozlamalari
st.set_page_config(layout="wide", page_title="Evapotranspiration Analysis")

# Google Earth Engine autentifikatsiya qilish
try:
    ee.Initialize()
except Exception as e:
    st.error("Google Earth Engine autentifikatsiya xatosi! Ilova qayta ishga tushirilishi kerak.")

# ---- 1. Header ----
st.title("🌍 Evapotranspiration Analysis")

# Sidebar sozlamalari
with st.sidebar:
    st.header("Sozlamalar")
    start_date = st.date_input("Boshlang'ich sana", datetime(2024, 1, 1))
    end_date = st.date_input("Tugash sana", datetime(2024, 2, 1))
    selected_result = st.selectbox("Natijani tanlang", ["NDVI", "ET", "Rn", "H", "G", "LAI"])

# ---- 2. Xarita yaratish ----
st.subheader("Interaktiv xarita")
Map = geemap.Map()
Map.add_basemap("HYBRID")

# Polygon chizish uchun funksiya
def get_polygon():
    feature = Map.user_roi  # Foydalanuvchi chizgan polygon
    if feature:
        return feature.geometry()
    else:
        st.warning("Iltimos, polygon chizing!")
        return None

Map.to_streamlit()

# ---- 3. Tugmani bosgandan keyin hisoblash ----
if st.button("Hisoblash"):
    polygon = get_polygon()
    if polygon:
        # Landsat 8 kolleksiyasini yuklash
        dataset = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA") \
            .filterBounds(polygon) \
            .filterDate(ee.Date(start_date.strftime('%Y-%m-%d')), ee.Date(end_date.strftime('%Y-%m-%d')))

        # Tanlangan natijaga qarab hisoblash
        if selected_result == "NDVI":
            def calculate_ndvi(image):
                ndvi = image.normalizedDifference(["B5", "B4"]).rename("NDVI")
                return image.addBands(ndvi)
            dataset = dataset.map(calculate_ndvi)
            final_image = dataset.median().select("NDVI")
            Map.addLayer(final_image, {"min": -1, "max": 1, "palette": ["blue", "white", "green"]}, "NDVI")

        # Xaritada natijani ko'rsatish
        st.success(f"{selected_result} natijasi xaritaga qo'shildi!")
        Map.to_streamlit()

# ---- 4. Natijani yuklab olish ----
if st.button("Natijani yuklab olish"):
    if selected_result == "NDVI":
        url = final_image.getDownloadURL()
        st.markdown(f'[Yuklab olish]( {url} )')



#def result_final(image,srtm,air_temp,wind_image,bosim):
    """Yakuniy natijani hisoblaymiz Azizlar !"""
    spek_rad = bk.calc_spek_rad(image) 
    reflictivity = bk.reflectivity(spek_rad)
    top_albedo = bk.albedo_toa(reflictivity)
    albedo = bk.calc_albedo(top_albedo)
    in_short = bk.calc_inshort_radiation(albedo)
    emissivity = bk.calc_surface_ems(in_short)
    surface_temp = bk.calc_surface_temp(emissivity)
    out_rad = bk.cal_out_rad(surface_temp)
    in_long = bk.calc_in_long(out_rad,air_temp, srtm)
    net_radian = bk.net_radiation(in_long)
    soil_heat = bk.soil_heat_flux(net_radian)
    aero_res = bk.calc_aerodinamic_res(soil_heat,wind_image)
    different_t = bk.different_temp(aero_res, air_temp)
    air_pressure = bk.calc_havo_bosim(different_t, bosim, wind_image)
    heat_flux = bk.calc_heat_flux(air_pressure)
    e_t = bk.evapotranspiration(heat_flux)
   
