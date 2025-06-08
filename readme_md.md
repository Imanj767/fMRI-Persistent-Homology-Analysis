# fMRI Persistent Homology Analysis App

این اپ برای محاسبه persistent homology از داده‌های fMRI در پلتفرم Brainlife طراحی شده است.

## توضیحات

این اپ مراحل زیر را انجام می‌دهد:

1. **بارگذاری داده‌های fMRI**: داده‌های fMRI پیش‌پردازش شده را از فرمت NIfTI خوانده می‌شود
2. **استخراج سری زمانی**: سری‌های زمانی از نواحی مغزی با استفاده از atlas یا روش خودکار استخراج می‌شود
3. **محاسبه اتصال عملکردی**: ماتریس اتصال عملکردی با روش‌های مختلف (correlation, covariance, precision) محاسبه می‌شود
4. **محاسبه Persistent Homology**: با استفاده از کتابخانه Ripser، persistent homology محاسبه می‌شود
5. **استخراج ویژگی**: ویژگی‌های مختلف از persistence diagrams استخراج می‌شود
6. **تجسم نتایج**: persistence diagrams و ماتریس اتصال تجسم می‌شود

## ورودی‌ها

### اجباری
- **fmri**: فایل داده‌های fMRI پیش‌پردازش شده (فرمت NIfTI)

### اختیاری
- **atlas**: فایل atlas برای تعریف نواحی مغزی (فرمت NIfTI)
- **n_regions**: تعداد نواحی برای تقسیم‌بندی خودکار (پیش‌فرض: 100)
- **maxdim**: حداکثر بعد برای محاسبه persistent homology (پیش‌فرض: 1)
- **connectivity_method**: روش محاسبه اتصال عملکردی (پیش‌فرض: "correlation")
- **low_pass**: فرکانس بالا گذر برای فیلتر (پیش‌فرض: 0.1)
- **high_pass**: فرکانس پایین گذر برای فیلتر (پیش‌فرض: 0.01)
- **t_r**: زمان تکرار (TR) به ثانیه (پیش‌فرض: 2.0)

## خروجی‌ها

### فایل‌های عددی
- `persistence_diagram_dim0.npy`: Persistence diagram برای بعد 0 (connected components)
- `persistence_diagram_dim1.npy`: Persistence diagram برای بعد 1 (loops)
- `persistence_features.json`: ویژگی‌های استخراج شده از persistence diagrams
- `persistence_features.csv`: ویژگی‌ها در فرمت CSV
- `connectivity_matrix.npy`: ماتریس اتصال عملکردی

### تصاویر
- `persistence_diagrams.png`: تجسم persistence diagrams
- `connectivity_matrix.png`: تجسم ماتریس اتصال عملکردی

### ویژگی‌های محاسبه شده

برای هر بعد persistent homology:
- `num_features`: تعداد ویژگی‌های topological  
- `max_persistence`: حداکثر persistence (طول عمر)
- `mean_persistence`: میانگین persistence
- `std_persistence`: انحراف معیار persistence
- `sum_persistence`: مجموع persistence
- `mean_birth`: میانگین زمان تولد
- `mean_death`: میانگین زمان مرگ
- `max_birth`: حداکثر زمان تولد
- `max_death`: حداکثر زمان مرگ

## نحوه استفاده

### در Brainlife
1. داده‌های fMRI پیش‌پردازش شده خود را آپلود کنید
2. اپ "fMRI Persistent Homology" را انتخاب کنید
3. پارامترهای مورد نظر را تنظیم کنید
4. اپ را اجرا کنید

### محلی (با Singularity)
```bash
# اجرای اپ (config.json باید موجود باشد)
./main

# یا مستقیماً:
singularity exec -e docker://brainlife/ga-python:lab328-dipy141-pybrainlife-1.0 python ./app.py
```

### محلی (بدون Singularity)
```bash
# نصب وابستگی‌ها
pip install -r requirements.txt

# اجرای اپ مستقیماً
python app.py
```

## مثال config.json

```json
{
  "fmri": "/path/to/preprocessed_fmri.nii.gz",
  "atlas": "/path/to/atlas.nii.gz",
  "n_regions": 100,
  "maxdim": 1,
  "connectivity_method": "correlation"
}
```

## الگوریتم Persistent Homology

Persistent homology یک ابزار قدرتمند در topological data analysis است که ساختارهای topological پایدار در داده‌ها را شناسایی می‌کند. در این اپ:

1. **ماتریس فاصله**: از ماتریس اتصال عملکردی، ماتریس فاصله محاسبه می‌شود
2. **Filtration**: یک مجموعه simplicial complexes با پارامتر فاصله ساخته می‌شود
3. **Homology**: برای هر مرحله، homology groups محاسبه می‌شود
4. **Persistence**: تغییرات homology در طول filtration رصد می‌شود

## کاربردهای علمی

- **تشخیص بیماری**: تشخیص تفاوت‌های topological بین افراد سالم و بیمار
- **تحلیل network**: درک ساختار پیچیده شبکه‌های مغزی
- **biomarker discovery**: کشف نشانگرهای زیستی جدید
- **developmental studies**: مطالعه تغییرات topological در طول رشد

## مراجع

1. Sporns, O. (2018). Graph theory methods: applications in brain networks. Dialogues in clinical neuroscience, 20(2), 111-121.
2. Petri, G., et al. (2014). Homological scaffolds of brain functional networks. Journal of the Royal Society Interface, 11(101), 20140873.
3. Chung, M. K., et al. (2018). Persistent homology in sparse regression and its application to brain morphometry. IEEE transactions on medical imaging, 34(9), 1928-1939.

## نویسنده

[نام شما]
[ایمیل شما]

## مجوز

MIT License
